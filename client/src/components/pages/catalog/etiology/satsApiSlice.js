import { catalogApiSlice } from '../../../../services/store/rootApi';
import SATSSignaturePresence from '../../../controls/plotly/SATS/satsSignaturePresence';
import SATSDotPlot from '../../../controls/plotly/SATS/satsDotPlot';
import { groupBy } from 'lodash';
import satsExampleData_SBS from '../../../controls/plotly/SATS/sats_example_data_SBS.json';

export const satsApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    satsDataBySignature: builder.query({
      query: (params) => ({
        url: 'signature_etiology',  // Use signature_etiology for SATS data
        params: { 
          signatureName: params.signatureName,
          signatureSetName: params.signatureSetName,
          limit: 1000000  // Get all data for comprehensive analysis
        },
      }),
      transformResponse: (data, meta, args) => {
        console.log('🎯 SATS API Raw Data:', {
          endpoint: 'signature_etiology',
          params: args,
          dataLength: data?.length || 0,
          sampleData: data?.slice(0, 3) || [],
          fields: data?.length > 0 ? Object.keys(data[0]) : []
        });

        if (!data || data.length === 0) {
          console.warn('⚠️ No data returned from signature_etiology API');
          return { traces: [], layout: {}, config: {} };
        }

        // Transform the etiology data into SATS format
        const transformedData = transformEtiologyDataToSATS(data);
        
        // Log the transformed data for debugging
        logSATSData(transformedData.tmbData);
        
        console.log('🎯 SATS Transformed Data:', {
          tmbDataLength: transformedData.tmbData?.length || 0,
          dotDataLength: transformedData.dotData?.length || 0,
          sampleTmbData: transformedData.tmbData?.slice(0, 3) || [],
          expectedFields: ['CancerType', 'SBS', 'Count', 'Proportion', 'N', 'Presence', 'TMB_all', 'CancerType_num']
        });
        
        // Generate the SATS plot
        return SATSSignaturePresence(transformedData);
      },
    }),
    satsExampleData: builder.query({
      queryFn: async () => {
        try {
          console.log('🎯 Creating Complete SATS Plot from example data...');
          console.log('📊 Example data sample:', satsExampleData_SBS?.slice(0, 3));
          
          if (!satsExampleData_SBS || satsExampleData_SBS.length === 0) {
            throw new Error('Example data not available');
          }

          // Transform the example data to the format expected by SATSSignaturePresence
          const transformedData = satsExampleData_SBS.map(item => {
            // Extract base signature name (e.g., "SBS1(Deamination of 5meC)" → "SBS1")
            const baseSignature = item.SBS.match(/^([^(]+)/)?.[1] || item.SBS;
            
            return {
              cancer: item.CancerType,
              signature: baseSignature.trim(), // Use extracted base signature
              presence: item.Presence, // Use Presence (capital P) from the JSON
              tmb: item.TMB_all || (item.Presence * 10), // Use TMB_all or scale Presence
              sampleCount: item.N || 100 // Use N (sample count) from data
            };
          });

          console.log('🔄 Transformed data sample:', transformedData.slice(0, 3));
          console.log('📊 Total records transformed:', transformedData.length);
          console.log('🏥 Unique cancer types:', [...new Set(transformedData.map(d => d.cancer))]);
          console.log('✍️ Unique signatures:', [...new Set(transformedData.map(d => d.signature))]);

          // Check if any records have presence > 0
          const validRecords = transformedData.filter(d => d.presence > 0);
          console.log('✅ Records with presence > 0:', validRecords.length);
          console.log('📈 Sample valid records:', validRecords.slice(0, 5));

          // Generate the complete SATS plot using SATSSignaturePresence function (TMB + Dot plot)
          console.log('🎯 Calling SATSSignaturePresence for complete plot...');
          const plotResult = SATSSignaturePresence({ 
            tmbData: transformedData, 
            dotData: transformedData 
          });
          
          console.log('🎯 SATSSignaturePresence result:', plotResult);
          console.log('📈 Plot traces:', plotResult?.traces);
          console.log('📊 Plot traces count:', plotResult?.traces?.length);
          console.log('🎨 Plot layout:', plotResult?.layout);
          
          // If SATSSignaturePresence doesn't work, fall back to just the dot plot
          if (!plotResult?.traces || plotResult.traces.length === 0) {
            console.log('⚠️ SATSSignaturePresence failed, falling back to dot plot only');
            const dotPlotResult = SATSDotPlot(transformedData);
            return { data: dotPlotResult };
          }
          
          return { data: plotResult };
        } catch (error) {
          console.error('❌ Error creating complete SATS plot:', error);
          
          // Final fallback: create simple dot plot
          try {
            console.log('🔄 Attempting fallback to dot plot...');
            const transformedData = satsExampleData.map(item => {
              const baseSignature = item.SBS.match(/^([^(]+)/)?.[1] || item.SBS;
              return {
                cancer: item.CancerType,
                signature: baseSignature.trim(),
                presence: item.Presence,
                tmb: item.TMB_all || (item.Presence * 10),
                sampleCount: item.N || 100
              };
            });
            
            const fallbackResult = SATSDotPlot(transformedData);
            return { data: fallbackResult };
          } catch (fallbackError) {
            console.error('❌ Fallback also failed:', fallbackError);
            return { error: { message: 'Failed to create SATS plot: ' + error.message } };
          }
        }
      },
    }),
    satsEtiologyLookup: builder.query({
      query: (params) => ({
        url: 'signature_etiology',
        params: { ...params, limit: 1000 },
      }),
      // Return raw data for signature lookup
      transformResponse: (data) => data,
    }),
  }),
});

// Transform function to convert etiology API data to SATS format
function transformEtiologyDataToSATS(etiologyData) {
  if (!etiologyData || etiologyData.length === 0) {
    return { tmbData: [], dotData: [] };
  }

  // Signature name to annotation mapping (from R code)
  const signatureAnnotations = {
    'SBS1': 'SBS1(Deamination of 5meC)',
    'SBS2': 'SBS2/13(APOBEC)',
    'SBS13': 'SBS2/13(APOBEC)',
    'SBS2_13': 'SBS2/13(APOBEC)',
    'SBS4': 'SBS4(Tobacco smoking)',
    'SBS6': 'SBS6(Defective MMR)',
    'SBS7a': 'SBS7a(UV exposure)',
    'SBS7b': 'SBS7b(UV exposure)',
    'SBS8': 'SBS8(Unknown)',
    'SBS10a': 'SBS10a(POLE-exo*)',
    'SBS10b': 'SBS10b(POLE-exo*)',
    'SBS10c': 'SBS10c(POLE-exo*)',
    'SBS11': 'SBS11(TMZ treatment)',
    'SBS12': 'SBS12(Unknown)',
    'SBS14': 'SBS14(Defective MMR)',
    'SBS15': 'SBS15(Defective MMR)',
    'SBS19': 'SBS19(Unknown)',
    'SBS25': 'SBS25(chemotherapy)',
    'SBS30': 'SBS30(Defective BER)',
    'SBS32': 'SBS32(AZA treatment)',
    'SBS44': 'SBS44(Defective MMR)',
    'SBS84': 'SBS84(AID)',
    'SBS87': 'SBS87(TP treatment)',
    'SBS89': 'SBS89(Unknown)',
    'SBS94': 'SBS94(Unknown)',
    'SBS97': 'SBS97(Unknown)',
    'SBS_Flat': 'Flat(SBS3/5/40a/40b)',
    'SBS_Artefactes': 'Artefactes(SBS50/51/57)'
  };

  // Function to get annotated signature name
  function getAnnotatedSignatureName(signatureName) {
    return signatureAnnotations[signatureName] || signatureName;
  }

  // Group data by cancer type first
  const groupedByCancer = groupBy(etiologyData, 'cancer');
  
  // Step 1: Calculate TMB data (equivalent to R TMB.all creation)
  const tmbDataMap = new Map();
  const cancerTMBTotals = new Map();
  
  Object.entries(groupedByCancer).forEach(([cancerType, cancerData]) => {
    // Group by signature for this cancer type
    const groupedBySignature = groupBy(cancerData, 'signatureName');
    
    // Get unique samples for this cancer type
    const uniqueSamples = [...new Set(cancerData.map(item => item.sample))];
    const totalSamples = uniqueSamples.length;
    
    let cancerTotalTMB = 0;
    
    Object.entries(groupedBySignature).forEach(([signatureName, signatureData]) => {
      // Calculate total exposure for this signature across all samples
      // This mimics the R code: sum(exp_sig_clinical) / assaySize
      const totalExposure = signatureData.reduce((sum, item) => sum + (item.exposure || 0), 0);
      
      // If burden field is available, use it; otherwise calculate from exposure
      const signatureTMB = signatureData.length > 0 && signatureData[0].burden 
        ? signatureData.reduce((sum, item) => sum + (item.burden || 0), 0) / signatureData.length
        : totalExposure / totalSamples; // Normalize by sample count
      
      cancerTotalTMB += signatureTMB;
      
      // Store TMB data
      const key = `${cancerType}_${signatureName}`;
      tmbDataMap.set(key, {
        CancerType: cancerType,
        SBS: signatureName,
        Count: signatureTMB, // TMB per signature (mutations per Mb)
        N: totalSamples,
        TMB_all: 0 // Will be filled in next step
      });
    });
    
    cancerTMBTotals.set(cancerType, cancerTotalTMB);
  });
  
  // Step 2: Update TMB_all for all records and sort cancer types by total TMB
  const sortedCancers = Array.from(cancerTMBTotals.entries())
    .sort((a, b) => b[1] - a[1]) // Sort by TMB descending (like R code)
    .map(([cancer, _]) => cancer);
    
  // Step 3: Calculate presence data and assign cancer type numbers
  const satsData = [];
  const cancerTypeMapping = new Map();
  
  sortedCancers.forEach((cancerType, index) => {
    cancerTypeMapping.set(cancerType, index + 1);
  });
  
  Object.entries(groupedByCancer).forEach(([cancerType, cancerData]) => {
    const groupedBySignature = groupBy(cancerData, 'signatureName');
    const uniqueSamples = [...new Set(cancerData.map(item => item.sample))];
    const totalSamples = uniqueSamples.length;
    const tmbAll = cancerTMBTotals.get(cancerType);
    const cancerTypeNum = cancerTypeMapping.get(cancerType);
    
    Object.entries(groupedBySignature).forEach(([signatureName, signatureData]) => {
      // Count samples with detectable signature (exposure > 0)
      // This mimics the R presence calculation
      const samplesWithSignature = signatureData.filter(item => (item.exposure || 0) > 0).length;
      const presence = totalSamples > 0 ? samplesWithSignature / totalSamples : 0;
      
      // Calculate proportion (signature's contribution to total TMB)
      const tmbKey = `${cancerType}_${signatureName}`;
      const tmbData = tmbDataMap.get(tmbKey);
      const proportion = tmbAll > 0 ? (tmbData?.Count || 0) / tmbAll : 0;
      
      satsData.push({
        CancerType: cancerType,
        SBS: getAnnotatedSignatureName(signatureName), // Use annotated signature name
        Count: tmbData?.Count || 0, // TMB per signature
        Proportion: proportion,
        N: totalSamples,
        Presence: presence,
        TMB_all: tmbAll,
        Label: totalSamples.toString(), // String version for display
        CancerType_num: cancerTypeNum
      });
    });
  });
  
  return {
    tmbData: satsData,
    dotData: satsData
  };
}

// Log the transformed data for debugging
function logSATSData(satsData) {
  if (satsData && satsData.length > 0) {
    console.log('🎯 SATS Data Sample (first 5 rows):');
    console.table(satsData.slice(0, 5));
    
    // Verify expected columns exist
    const expectedColumns = ['CancerType', 'SBS', 'Count', 'Proportion', 'N', 'Presence', 'TMB_all', 'Label', 'CancerType_num'];
    const actualColumns = Object.keys(satsData[0] || {});
    const missingColumns = expectedColumns.filter(col => !actualColumns.includes(col));
    
    if (missingColumns.length > 0) {
      console.warn('⚠️ Missing expected columns:', missingColumns);
    } else {
      console.log('✅ All expected columns present:', expectedColumns);
    }
  }
}

// Alternative transformation if you have different data structure
function transformSignatureActivityToSATS(activityData) {
  // This would be used if you're getting data from signature_activity endpoint
  const groupedData = groupBy(activityData, (item) => 
    `${item.signatureName}_${item.cancer}`
  );

  const tmbData = [];
  const dotData = [];

  Object.values(groupedData).forEach(group => {
    const firstItem = group[0];
    const { signatureName, cancer } = firstItem;

    // For TMB: use cancerBurden if available, or calculate from exposure
    const tmbValues = group
      .filter(item => item.cancerBurden && item.cancerBurden > 0)
      .map(item => item.cancerBurden);
    
    if (tmbValues.length > 0) {
      const avgTMB = tmbValues.reduce((sum, val) => sum + val, 0) / tmbValues.length;
      
      tmbData.push({
        cancer,
        signature: signatureName,
        tmb: avgTMB,
        sampleCount: group.length
      });
    }

    // For presence: proportion of samples with exposure > threshold
    const threshold = 0.01; // 1% exposure threshold
    const totalSamples = group.length;
    const samplesWithSignature = group.filter(item => 
      item.exposure && item.exposure > threshold
    ).length;
    
    const presence = totalSamples > 0 ? samplesWithSignature / totalSamples : 0;
    
    if (presence > 0) {
      dotData.push({
        cancer,
        signature: signatureName,
        presence,
        sampleCount: totalSamples
      });
    }
  });

  return {
    tmbData,
    dotData
  };
}

export const { useSatsDataBySignatureQuery, useSatsExampleDataQuery, useSatsEtiologyLookupQuery } = satsApiSlice;

// Export transformation functions for reuse
export { transformEtiologyDataToSATS, transformSignatureActivityToSATS };
