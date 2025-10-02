import React from 'react';

export default function Instructions() {
  return (
    <div className="bg-white border rounded p-3">
      <h4>Refitting Instructions</h4>
      <p>
        Welcome to the Signature Refitting module. This tool allows you to 
        comprehensively evaluate the accuracy of mutational signatures based on 
        different statistical variables (Cosine similarity, BIC, L2 norms etc) 
        and re-decompose signatures using different algorithms (SigProfiler, 
        deconstructsig, bootstrapping method).
      </p>
      
      <h5>How to Use:</h5>
      <ol>
        <li>Upload your mutation data or select from example datasets</li>
        <li>Choose the reference signatures to refit against</li>
        <li>Select the statistical metrics you want to evaluate</li>
        <li>Configure the refitting parameters</li>
        <li>Submit your analysis</li>
        <li>Monitor progress in the Status tab</li>
        <li>View results in the Targeted Sequencing tab</li>
      </ol>

      <h5>Supported File Formats:</h5>
      <ul>
        <li>VCF files</li>
        <li>MAF files</li>
        <li>TSV/CSV mutation matrices</li>
      </ul>

      <h5>Statistical Metrics:</h5>
      <ul>
        <li><strong>Cosine Similarity:</strong> Measures the similarity between original and refitted signatures</li>
        <li><strong>BIC (Bayesian Information Criterion):</strong> Model selection criterion</li>
        <li><strong>L2 Norms:</strong> Euclidean distance measure</li>
        <li><strong>Reconstruction Error:</strong> Difference between original and reconstructed data</li>
      </ul>
    </div>
  );
}
