const express = require('express');
const { v4: uuidv4 } = require('uuid');
const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');
const os = require('os');
const { del } = require('@vercel/blob');
const { pipeline } = require('node:stream/promises');
const pLimit = require('p-limit');
require('dotenv').config();

const app = express();
const port = process.env.PORT || 8000;
app.use(express.json());

// ====================================================================
// HELPER FUNCTIONS
// ====================================================================

/**
 * Runs an R script as a synchronous child process.
 * @param {string} scriptName - The name of the R script in the /scripts folder.
 * @param {string} inputFilePath - The path (or comma-separated paths) to the input data file(s).
 * @returns {object | undefined} - The parsed JSON object from the R script's output file, or undefined if it fails.
 */
function runRScriptSync(scriptName, inputFilePath) {
    const scriptPath = path.join(__dirname, 'scripts', scriptName);
    const outputFilePath = path.join(os.tmpdir(), `${uuidv4()}.json`);
    
    // inputFilePath might be a comma-separated string of paths.
    // We wrap it in quotes to handle spaces, but R script logic parses it.
    const command = `Rscript "${scriptPath}" "${inputFilePath}" "${outputFilePath}"`;
    
    console.log(`Executing command: ${command}`);
    try {
        execSync(command, { stdio: 'pipe', timeout: 600000 }); // 10 minute timeout
        if (fs.existsSync(outputFilePath)) {
            const resultJson = fs.readFileSync(outputFilePath, 'utf-8');
            fs.unlinkSync(outputFilePath);
            return JSON.parse(resultJson);
        } else {
            console.error(`R script '${scriptName}' finished but did not create an output file.`);
            return undefined;
        }
    } catch (error) {
        const stderr = error.stderr ? error.stderr.toString() : 'No stderr output.';
        console.error('R SCRIPT STDERR:', stderr);
        if (fs.existsSync(outputFilePath)) {
            fs.unlinkSync(outputFilePath);
        }
        throw new Error(stderr || `R script '${scriptName}' failed.`);
    }
}

/**
 * Downloads a file from a URL to a local path using native fetch.
 * @param {string} fileUrl - The public URL of the file.
 * @param {string} outputPath - The local path to save the file.
 */
async function downloadFile(fileUrl, outputPath) {
    const response = await fetch(fileUrl);
    if (!response.ok) {
        throw new Error(`Failed to download file from ${fileUrl}: ${response.status} ${response.statusText}`);
    }
    await pipeline(response.body, fs.createWriteStream(outputPath));
}

/**
 * Deletes files from Vercel Blob storage.
 * @param {string|string[]} fileUrls - A single URL or array/comma-separated string of URLs.
 */
async function deleteFromBlob(fileUrls) {
    // Convert comma-separated string or array to array
    const urlsToDelete = Array.isArray(fileUrls) ? fileUrls : fileUrls.split(',');
    
    for (const url of urlsToDelete) {
        if (!url) continue;
        try {
            await del(url.trim(), { token: process.env.BLOB_READ_WRITE_TOKEN });
            console.log(`Successfully deleted from Blob: ${url}`);
        } catch (error) {
            console.error(`Failed to delete from Blob: ${url}`, error);
        }
    }
}

/**
 * Sends the final job status back to the Vercel API.
 */
async function updateVercelStatus(jobId, status, data) {
    console.log(`[${jobId}] Updating Vercel status to: ${status}`);
    try {
        const response = await fetch(process.env.VERCEL_STATUS_UPDATE_URL, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ jobId, status, ...data }),
        });
        if (!response.ok) {
          console.error(`[${jobId}] Vercel status update failed with status: ${response.status}`);
        }
    } catch (err) {
        console.error(`[${jobId}] CRITICAL: FAILED to update Vercel status. Error:`, err.message);
    }
}

/**
 * Identifies a single GC-MS peak using MoNA.
 */
async function identifySinglePeak(peakData) {
    // 1. Setup a timeout controller (fails request after 10 seconds)
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 10000);

    try {
        const spectrumString = peakData.spectrum_data
            .map(p => `${p.mz.toFixed(4)}:${(p.relative_intensity || 0).toFixed(0)}`)
            .join(" ");
            
        const payload = { spectrum: spectrumString, minSimilarity: 500, algorithm: "default" };
        
        const response = await fetch("https://mona.fiehnlab.ucdavis.edu/rest/similarity/search", {
            method: 'POST',
            headers: { 
                'Content-Type': 'application/json',
                'User-Agent': 'MassSpecWorker/1.0' // Good practice to avoid blocks
            },
            body: JSON.stringify(payload),
            signal: controller.signal // 2. Attach the abort signal here
        });

        if (!response.ok) {
            // If API returns 500, 404, etc., just return null
            return null;
        }

        const apiResults = await response.json();
        
        if (apiResults && apiResults.length > 0) {
            const bestHit = apiResults[0];
            const compoundName = bestHit.hit?.compound?.[0]?.names?.[0]?.name || "Unknown";
            const librarySpectrum = bestHit.hit?.spectrum ? bestHit.hit.spectrum.split(' ').map(p => {
                const parts = p.split(':');
                return { mz: parseFloat(parts[0]), intensity: parseFloat(parts[1]) };
            }) : [];

            return { 
                peak_number: peakData.peak_number, 
                match_name: compoundName, 
                similarity_score: bestHit.score,
                library_spectrum: librarySpectrum
            };
        }
    } catch (error) {
        // 3. Log the error simply, and return null so the job doesn't crash
        console.warn(`Peak #${peakData.peak_number} failed/timeout: ${error.message}`);
        return null; 
    } finally {
        // 4. Always clear the timer to prevent memory leaks
        clearTimeout(timeoutId);
    }

    // Return null if no matches found
    return null;
}


// ====================================================================
// MAIN WORKER ENDPOINT
// ====================================================================
app.post('/process-job', async (req, res) => {
    if (req.headers['x-worker-secret'] !== process.env.WORKER_SECRET) {
        return res.status(401).send('Unauthorized');
    }

    // `fileUrl` can now be a comma-separated string (e.g. "url1,url2") from the updated UI
    const { jobId, analysisType, fileUrl, originalFilename, inputDataHash } = req.body;
    
    console.log(`[${jobId}] Received job. Type: ${analysisType}`);
    res.status(202).send('Accepted');

    // Parse inputs (handles single or multiple files)
    const fileUrls = fileUrl.split(',');
    const fileNames = originalFilename.split(',');
    const downloadedFilePaths = [];

    try {
        // --- DOWNLOAD STAGE ---
        for (let i = 0; i < fileUrls.length; i++) {
            const url = fileUrls[i].trim();
            const name = fileNames[i] ? fileNames[i].trim() : `file_${i}.mzML`;
            const localPath = path.join(os.tmpdir(), `${uuidv4()}_${name}`);
            
            console.log(`[${jobId}] Downloading file ${i + 1}/${fileUrls.length}: ${name}`);
            await downloadFile(url, localPath);
            downloadedFilePaths.push(localPath);
        }

        const scriptMap = { 'xcms': 'xcms_analysis.R', 'drc': 'drc_analysis.R', 'nmr': 'nmr1d_analysis.R' };
        const scriptName = scriptMap[analysisType];
        if (!scriptName) throw new Error(`Invalid analysis type: ${analysisType}`);

        // Construct input argument for R (comma-separated local paths)
        const rScriptInput = downloadedFilePaths.join(',');

        // --- STAGE 1: Run R Script ---
        console.log(`[${jobId}] Running R script on ${downloadedFilePaths.length} file(s)...`);
        const rScriptResults = runRScriptSync(scriptName, rScriptInput);

        if (!rScriptResults || typeof rScriptResults !== 'object') {
            throw new Error('R script finished but did not produce a valid result object.');
        }
        
        if (rScriptResults.status !== 'success') {
            throw new Error(rScriptResults.error || 'The R script reported a failure.');
        }
        console.log(`[${jobId}] Stage 1 (R Script) completed successfully.`);
        
        // --- STAGE 2: Run MoNA Identification (only for 'xcms' jobs) ---
        if (analysisType === 'xcms') {
            console.log(`[${jobId}] Stage 2 (MoNA Identification) starting...`);
            
            // 1. Get the full list from R results
            const rawSpectra = rScriptResults.results?.top_spectra_data || [];
            
            // 2. FORCE LIMIT: Explicitly take only the first 50 items
            // This protects your server if the R script returns 1000 peaks
            const topSpectra = rawSpectra.slice(0, 50); 

            console.log(`[${jobId}] Processing ${topSpectra.length} spectra (filtered from ${rawSpectra.length} total)...`);

            const limit = pLimit.default(5);
            
            if (topSpectra.length > 0) {
                const allMatches = await Promise.all(
                    topSpectra.map(spec => limit(() => identifySinglePeak(spec)))
                );
                rScriptResults.results.library_matches = allMatches.filter(Boolean);
            } else {
                rScriptResults.results.library_matches = [];
            }
            console.log(`[${jobId}] Stage 2 (MoNA Identification) complete.`);
        }
        
        await updateVercelStatus(jobId, 'completed', { result: rScriptResults });

    } catch (error) {
        console.error(`[${jobId}] JOB FAILED:`, error.message);
        await updateVercelStatus(jobId, 'failed', { error: error.message });
    
    } finally {
        // Cleanup local temp files
        for (const filePath of downloadedFilePaths) {
            if (fs.existsSync(filePath)) {
                fs.unlinkSync(filePath);
            }
        }
        console.log(`[${jobId}] Cleaned up local temp files.`);
        
        // Cleanup Blob storage
        await deleteFromBlob(fileUrl);
    }
});

app.get('/healthcheck', (req, res) => res.json({ status: "OK" }));

app.listen(port, () => {
    console.log(`Asynchronous Node.js Worker running on port ${port}`);
});