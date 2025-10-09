// server.js - Synchronous Node.js API for R Scripts (Updated Paths)

const express = require('express');
const fileUpload = require('express-fileupload');
const cors = require('cors');
const { v4: uuidv4 } = require('uuid');
const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');
const os = require('os');

const app = express();
const port = 8000;
const MAX_FILE_SIZE = 10 * 1024 * 1024;

app.use(cors());
app.use(fileUpload({ limits: { fileSize: MAX_FILE_SIZE } }));
app.use(express.text({ type: ['text/csv', 'text/plain'] }));
app.use(express.json())
// --- THIS IS YOUR ORIGINAL, UNTOUCHED HELPER FUNCTION ---
function runRScriptSync(scriptName, inputFilePath) {
    const scriptPath = path.join(__dirname, 'scripts', scriptName);
    const outputFilePath = path.join(os.tmpdir(), `${uuidv4()}.json`);
    const command = `Rscript "${scriptPath}" "${inputFilePath}" "${outputFilePath}"`;
    console.log(`Executing command: ${command}`);

    try {
        execSync(command, { stdio: 'pipe', timeout: 300000 });
        const resultJson = fs.readFileSync(outputFilePath, 'utf-8');
        fs.unlinkSync(outputFilePath);
        return JSON.parse(resultJson);
    } catch (error) {
        console.error(`Error executing R script: ${scriptName}`);
        const stderr = error.stderr ? error.stderr.toString() : 'No stderr output.';
        console.error('STDERR:', stderr);
        if (fs.existsSync(outputFilePath)) {
            fs.unlinkSync(outputFilePath);
        }
        throw new Error(stderr);
    }
}

// --- START: NEW HELPER FUNCTION (ADDITION ONLY) ---
// This new function is specifically for the identifier script, which needs 3 arguments.
function runIdentifierScriptSync(scriptName, initialJsonPath, outputPath) {
    const scriptPath = path.join(__dirname, 'scripts', scriptName);
    // Note the three arguments being passed
    const command = `Rscript "${scriptPath}" "${initialJsonPath}" "${outputPath}"`;
    console.log(`Executing identifier command: ${command}`);

    try {
        execSync(command, { stdio: 'pipe', timeout: 300000 });
        const resultJson = fs.readFileSync(outputPath, 'utf-8');
        // No need to clean up here, the 'finally' block in the endpoint will do it.
        return JSON.parse(resultJson);
    } catch (error) {
        console.error(`Error executing R script: ${scriptName}`);
        const stderr = error.stderr ? error.stderr.toString() : 'No stderr output.';
        console.error('STDERR:', stderr);
        throw new Error(stderr);
    }
}
// --- END: NEW HELPER FUNCTION ---


app.get('/healthcheck', (req, res) => {
    res.json({ status: "OK", timestamp: new Date() });
});

// --- THIS IS YOUR ORIGINAL, UNTOUCHED XCMS ENDPOINT ---
app.post('/analyze/xcms', (req, res) => {
    if (!req.files || !req.files.file) {
        return res.status(400).json({ error: "No file was uploaded in a form field named 'file'." });
    }
    const file = req.files.file;
    
    const tempFilePath = path.join(os.tmpdir(), `${uuidv4()}_${file.name}`);
    
    file.mv(tempFilePath, (err) => {
        if (err) return res.status(500).json({ error: 'Failed to save uploaded file.', details: err });
        try {
            // It calls the original, simple helper function
            const results = runRScriptSync('xcms_analysis.R', tempFilePath); // Ensure script name is correct
            res.json(results);
        } catch (error) {
            res.status(500).json({ status: 'error', error: error.message });
        } finally {
            fs.unlinkSync(tempFilePath);
        }
    });
});


// --- START: NEW IDENTIFICATION ENDPOINT (ADDITION ONLY) ---
app.post('/analyze/xcms/identify', (req, res) => {
    // It no longer looks for a file, only the JSON body.
    if (!req.body || !req.body.initialAnalysisData) {
        return res.status(400).json({ error: "Missing 'initialAnalysisData' in the request body." });
    }
    const initialDataString = JSON.stringify(req.body.initialAnalysisData);
    
    const tempId = uuidv4();
    const tempInitialJsonPath = path.join(os.tmpdir(), `${tempId}_initial.json`);
    const tempOutputPath = path.join(os.tmpdir(), `${tempId}_matches.json`);

    try {
        // Write the data from the UI to a temporary file for the R script to read
        fs.writeFileSync(tempInitialJsonPath, initialDataString);
        
        // Call the new, simpler helper function
        const results = runIdentifierScriptSync('ms_identifier.R', tempInitialJsonPath, tempOutputPath);
        res.json(results);
    } catch (error) {
        res.status(500).json({ status: 'error', error: error.message });
    } finally {
        // Clean up both temporary files
        if (fs.existsSync(tempInitialJsonPath)) fs.unlinkSync(tempInitialJsonPath);
        if (fs.existsSync(tempOutputPath)) fs.unlinkSync(tempOutputPath);
    }
});
// --- END: NEW IDENTIFICATION ENDPOINT ---

app.post('/analyze/drc', (req, res) => {
    const rawCsvData = req.body;
    if (!rawCsvData || typeof rawCsvData !== 'string' || rawCsvData.length === 0) {
        return res.status(400).json({ error: "Request body must contain raw CSV data." });
    }
    const tempFilePath = path.join(os.tmpdir(), `${uuidv4()}.csv`);
    try {
        fs.writeFileSync(tempFilePath, rawCsvData);
        const results = runRScriptSync('drc_analysis.R', tempFilePath);
        res.json(results);
    } catch (error) {
        res.status(500).json({ status: 'error', error: error.message });
    } finally {
        if (fs.existsSync(tempFilePath)) fs.unlinkSync(tempFilePath);
    }
});

app.listen(port, () => {
    console.log(`Synchronous Node.js API server running on http://localhost:${port}`);
});