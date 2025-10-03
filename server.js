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

app.use(cors());
app.use(fileUpload());
app.use(express.text({ type: ['text/csv', 'text/plain'] }));

function runRScriptSync(scriptName, inputFilePath) {
    const scriptPath = path.join(__dirname, 'scripts', scriptName);
    // 1. Create a unique path for the R script's output file
    const outputFilePath = path.join(os.tmpdir(), `${uuidv4()}.json`);

    // 2. Pass BOTH paths to the R script
    const command = `Rscript "${scriptPath}" "${inputFilePath}" "${outputFilePath}"`;

    console.log(`Executing command: ${command}`);

    try {
        // 3. Execute the command. We no longer need to capture stdout.
        // We set a very large timeout (5 minutes) for long analyses.
        execSync(command, { stdio: 'pipe', timeout: 300000 });

        // 4. Read the result from the file that the R script created.
        const resultJson = fs.readFileSync(outputFilePath, 'utf-8');
        
        // 5. Clean up the output file.
        fs.unlinkSync(outputFilePath);

        return JSON.parse(resultJson);
    } catch (error) {
        console.error(`Error executing R script: ${scriptName}`);
        console.error('STDERR:', error.stderr ? error.stderr.toString() : 'No stderr output.');
        
        // Clean up the output file even on error, if it exists
        if (fs.existsSync(outputFilePath)) {
            fs.unlinkSync(outputFilePath);
        }

        throw new Error(error.stderr ? error.stderr.toString() : 'The R script failed with an unknown error.');
    }
}

app.get('/healthcheck', (req, res) => {
    res.json({ status: "OK", timestamp: new Date(), wd: process.cwd() });
});

// --- Analysis Endpoints ---

app.post('/analyze/nmr', (req, res) => {
    if (!req.files || !req.files.file) {
        return res.status(400).json({ error: "No file was uploaded in a form field named 'file'." });
    }
    const file = req.files.file;
    const tempFilePath = path.join(os.tmpdir(), `${uuidv4()}.zip`);
    file.mv(tempFilePath, (err) => {
        if (err) return res.status(500).json({ error: 'Failed to save uploaded file.', details: err });
        try {
            console.log("calling")
            const results = runRScriptSync('nmr1d_analysis.R', tempFilePath);
            res.json(results);
        } catch (error) {
            console.log(error.message)
            res.status(500).json({ status: 'error', error: error.message });
        } finally {
            fs.unlinkSync(tempFilePath);
        }
    });
});

app.post('/analyze/xcms', (req, res) => {
    if (!req.files || !req.files.file) {
        return res.status(400).json({ error: "No file was uploaded in a form field named 'file'." });
    }
    const file = req.files.file;
    const tempFilePath = path.join(os.tmpdir(), `${uuidv4()}.mzML`);
    file.mv(tempFilePath, (err) => {
        if (err) return res.status(500).json({ error: 'Failed to save uploaded file.', details: err });
        try {
            const results = runRScriptSync('xcms_analysis.R', tempFilePath);
            res.json(results);
        } catch (error) {
            res.status(500).json({ status: 'error', error: error.message });
        } finally {
            fs.unlinkSync(tempFilePath);
        }
    });
});

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