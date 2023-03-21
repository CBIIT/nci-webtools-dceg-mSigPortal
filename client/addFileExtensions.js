// script for adding file extensions to local file imports
// may need adjustment for erroneously modified paths
const fs = require('fs/promises');
const path = require('path');
async function processFolder(dir, processFile) {
  const entries = await fs.readdir(dir, { withFileTypes: true });
  for (const entry of entries) {
    const fullPath = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      await processFolder(fullPath, processFile);
    } else if (entry.isFile() && path.extname(entry.name) === '.js') {
      await processFile(fullPath);
    }
  }
}
async function processFile(filePath) {
  const fileContent = await fs.readFile(filePath, 'utf-8');
  // Regex pattern to match relative ESM imports
  const importRegex = /import\s+(?:[^'"]+\s+from\s+)?['"](\.\/[^.'"]+)['"]/g;
  const updatedContent = fileContent.replace(
    importRegex,
    (match, importPath) => {
      // Add the '.js' extension to the import path
      const updatedImportPath = `${importPath}.js`;
      // Replace the original import path with the updated one
      return match.replace(importPath, updatedImportPath);
    }
  );
  // Write the updated content to the file
  await fs.writeFile(filePath, updatedContent, 'utf-8');
}
// Start the script by providing the folder path
const folderPath = './src/components/controls/plotly';
processFolder(folderPath, processFile)
  .then(() => console.log('File extensions added successfully.'))
  .catch((error) => console.error('An error occurred:', error));
