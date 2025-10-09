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


app.get('/healthcheck', (req, res) => {
    res.json({ status: "OK", timestamp: new Date() });
});

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

app.listen(port, () => {
    console.log(`Synchronous Node.js API server running on http://localhost:${port}`);
});