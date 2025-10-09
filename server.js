const express = require('express');
const { v4: uuidv4 } = require('uuid');
const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');
const os = require('os');
const { del } = require('@vercel/blob');
const { pipeline } = require('node:stream/promises');
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
 * @param {string} inputFilePath - The path to the input data file.
 * @returns {object | undefined} - The parsed JSON object from the R script's output file, or undefined if it fails.
 */
function runRScriptSync(scriptName, inputFilePath) {
    const scriptPath = path.join(__dirname, 'scripts', scriptName);
    const outputFilePath = path.join(os.tmpdir(), `${uuidv4()}.json`);
    const command = `Rscript "${scriptPath}" "${inputFilePath}" "${outputFilePath}"`;
    console.log(`Executing command: ${command}`);
    try {
        execSync(command, { stdio: 'pipe', timeout: 600000 }); // 10 minute timeout
        if (fs.existsSync(outputFilePath)) {
            const resultJson = fs.readFileSync(outputFilePath, 'utf-8');
            fs.unlinkSync(outputFilePath);
            return JSON.parse(resultJson);
        } else {
            // This case handles when the R script exits successfully but creates no output file.
            console.error(`R script '${scriptName}' finished but did not create an output file.`);
            return undefined;
        }
    } catch (error) {
        const stderr = error.stderr ? error.stderr.toString() : 'No stderr output.';
        console.error('R SCRIPT STDERR:', stderr);
        if (fs.existsSync(outputFilePath)) {
            fs.unlinkSync(outputFilePath);
        }
        // Re-throw the error to be caught by the main handler
        throw new Error(stderr || `R script '${scriptName}' failed.`);
    }
}

/**
 * Downloads a file from a URL to a local path using native fetch.
 * @param {string} fileUrl - The public URL of the file to download.
 * @param {string} outputPath - The local path to save the file to.
 */
async function downloadFile(fileUrl, outputPath) {
    const response = await fetch(fileUrl);
    if (!response.ok) {
        throw new Error(`Failed to download file: ${response.status} ${response.statusText}`);
    }
    await pipeline(response.body, fs.createWriteStream(outputPath));
}

/**
 * Deletes a file from Vercel Blob storage.
 * @param {string} fileUrl - The URL of the blob to delete.
 */
async function deleteFromBlob(fileUrl) {
    try {
        await del(fileUrl, { token: process.env.BLOB_READ_WRITE_TOKEN });
        console.log(`Successfully deleted from Blob: ${fileUrl}`);
    } catch (error) {
        console.error(`Failed to delete from Blob: ${fileUrl}`, error);
    }
}

/**
 * Sends the final job status back to the Vercel API.
 * @param {string} jobId - The ID of the job.
 * @param {'completed' | 'failed'} status - The final status.
 * @param {object} data - An object containing either a `result` or `error` property.
 */
async function updateVercelStatus(jobId, status, data) {
    console.log(`[${jobId}] Updating Vercel status to: ${status}`);
    try {
        const response = await fetch(process.env.VERCEL_STATUS_UPDATE_URL, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
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
 * Identifies a single GC-MS peak using the MoNA online database API.
 * @param {object} peakData - The peak data object from the XCMS script.
 * @returns {Promise<object>} - An object with the identification results.
 */
async function identifySinglePeak(peakData) {
    try {
        const spectrumString = peakData.spectrum_data.map(p => `${p.mz.toFixed(4)}:${(p.relative_intensity || 0).toFixed(0)}`).join(" ");
        const payload = { spectrum: spectrumString, minSimilarity: 500, algorithm: "default" };
        
        const response = await fetch("https://mona.fiehnlab.ucdavis.edu/rest/similarity/search", {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload),
        });

        if (!response.ok) return null;
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
        console.error(`Failed to identify peak #${peakData.peak_number}:`, error);
    }
    return { peak_number: peakData.peak_number, match_name: "No Match Found", similarity_score: 0, library_spectrum: [] };
}


// ====================================================================
// MAIN WORKER ENDPOINT
// ====================================================================
app.post('/process-job', async (req, res) => {
    if (req.headers['x-worker-secret'] !== process.env.WORKER_SECRET) {
        return res.status(401).send('Unauthorized');
    }

    const { jobId, analysisType, fileUrl, originalFilename, inputDataHash } = req.body;
    console.log(`[${jobId}] Received job. Type: ${analysisType}`);
    res.status(202).send('Accepted');

    const tempFilePath = path.join(os.tmpdir(), `${uuidv4()}_${originalFilename}`);

    try {
        await downloadFile(fileUrl, tempFilePath);
        console.log(`[${jobId}] File downloaded.`);

        const scriptMap = { 'xcms': 'xcms_analysis.R', 'drc': 'drc_analysis.R', 'nmr': 'nmr1d_analysis.R' };
        const scriptName = scriptMap[analysisType];
        if (!scriptName) throw new Error(`Invalid analysis type: ${analysisType}`);

        // --- STAGE 1: Run R Script ---
        const rScriptResults = runRScriptSync(scriptName, tempFilePath);

        // Add a robustness check to ensure the R script produced a valid result object.
        if (!rScriptResults || typeof rScriptResults !== 'object') {
            throw new Error('R script finished but did not produce a valid result object.');
        }
        
        // Safely check the status property.
        if (rScriptResults.status !== 'success') {
            throw new Error(rScriptResults.error || 'The R script reported a failure.');
        }
        console.log(`[${jobId}] Stage 1 (R Script) completed successfully.`);
        
        // --- STAGE 2: Run MoNA Identification (only for 'xcms' jobs) ---
        if (analysisType === 'xcms') {
            console.log(`[${jobId}] Stage 2 (MoNA Identification) starting...`);
            const topSpectra = rScriptResults.results?.top_spectra_data || [];
            
            if (topSpectra.length > 0) {
                const allMatches = await Promise.all(
                    topSpectra.map(spec => identifySinglePeak(spec))
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
        return;
    
    } finally {
        if (fs.existsSync(tempFilePath)) {
            fs.unlinkSync(tempFilePath);
            console.log(`[${jobId}] Cleaned up local temp file.`);
        }
    }

    await deleteFromBlob(fileUrl);
});

app.get('/healthcheck', (req, res) => res.json({ status: "OK" }));

app.listen(port, () => {
    console.log(`Asynchronous Node.js Worker running on port ${port}`);
});