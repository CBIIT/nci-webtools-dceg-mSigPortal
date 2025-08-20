// Browser Console Testing Script for SATS Plot
// Instructions: 
// 1. Open mSigPortal in your browser
// 2. Open Developer Tools (F12) 
// 3. Navigate to Catalog > Signature Catalog page
// 4. Paste this script in the console and run: testSATSIntegration()

window.testSATSIntegration = async function() {
  console.log('🧪 Testing SATS Integration...\n');
  
  try {
    console.log('1. Testing API endpoints...');
    
    // Test etiology options endpoint
    const optionsResponse = await fetch('/api/signature_etiology_options');
    const options = await optionsResponse.json();
    
    console.log('✅ Etiology options loaded:', {
      totalOptions: options.length,
      firstOption: options[0],
      availableStudies: [...new Set(options.map(o => o.study))].filter(Boolean),
      availableSignatureSets: [...new Set(options.map(o => o.signatureSetName))].filter(Boolean).slice(0, 5)
    });
    
    // Test etiology data endpoint with a sample query
    if (options.length > 0) {
      const sampleParams = new URLSearchParams({
        study: options[0].study || 'TCGA',
        signatureSetName: options[0].signatureSetName || 'COSMIC_v3.3_SBS_GRCh37',
        limit: 100
      });
      
      console.log('\n2. Testing etiology data endpoint...');
      console.log('🔗 Sample query:', `/api/signature_etiology?${sampleParams}`);
      
      const dataResponse = await fetch(`/api/signature_etiology?${sampleParams}`);
      const data = await dataResponse.json();
      
      console.log('✅ Etiology data loaded:', {
        totalRecords: data.length,
        sampleRecord: data[0],
        dataFields: data[0] ? Object.keys(data[0]) : [],
        uniqueCancers: [...new Set(data.map(d => d.cancer))].filter(Boolean).slice(0, 5),
        uniqueSignatures: [...new Set(data.map(d => d.signatureName))].filter(Boolean).slice(0, 5),
        burdenRange: data.length > 0 ? {
          min: Math.min(...data.map(d => d.burden || 0)),
          max: Math.max(...data.map(d => d.burden || 0)),
          average: (data.reduce((sum, d) => sum + (d.burden || 0), 0) / data.length).toFixed(2)
        } : 'No data'
      });
      
      // Test data transformation for SATS
      console.log('\n3. Testing SATS data transformation...');
      
      // Group by cancer type to see data structure
      const groupedByCancer = data.reduce((acc, curr) => {
        if (!acc[curr.cancer]) acc[curr.cancer] = [];
        acc[curr.cancer].push(curr);
        return acc;
      }, {});
      
      console.log('📊 Data grouped by cancer type:', {
        cancerTypes: Object.keys(groupedByCancer),
        sampleCancerData: Object.keys(groupedByCancer)[0] ? {
          cancer: Object.keys(groupedByCancer)[0],
          sampleCount: groupedByCancer[Object.keys(groupedByCancer)[0]].length,
          sampleRecords: groupedByCancer[Object.keys(groupedByCancer)[0]].slice(0, 3)
        } : 'No data'
      });
      
    } else {
      console.log('❌ No options data available to test with');
    }
    
    console.log('\n4. Navigation instructions:');
    console.log('🎯 To see the SATS plot:');
    console.log('   1. Navigate to: Catalog > Signature Catalog');
    console.log('   2. Select any category (e.g., "Cosmic Mutational Signatures")');
    console.log('   3. Scroll down to the "SATS - Signature Activity in Tumor Samples" section');
    console.log('   4. Select a Study and Signature Set');
    console.log('   5. Click "Generate SATS Plot"');
    
    console.log('\n✅ SATS integration test completed!');
    
  } catch (error) {
    console.error('❌ Error testing SATS integration:', error);
  }
};

// Auto-run the test
console.log('SATS Integration Tester loaded! Run: testSATSIntegration()');
