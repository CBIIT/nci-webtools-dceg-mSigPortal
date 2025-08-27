import { catalogApiSlice } from '../../../../services/store/rootApi';
import SATSSignaturePresence from '../../../controls/plotly/SATS/satsSignaturePresence';
import { groupBy } from 'lodash';

export const satsApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    satsOptions: builder.query({
      query: () => ({
        url: 'signature_etiology',
        params: { limit: 1000 },
      }),
      transformResponse: (data) => {
        if (!data || data.length === 0) {
          return [];
        }
        
        // Extract unique combinations of study and signatureSetName
        const uniqueCombinations = [];
        const seen = new Set();
        
        data.forEach(item => {
          const key = `${item.study}-${item.signatureSetName}`;
          if (!seen.has(key) && item.study && item.signatureSetName) {
            seen.add(key);
            uniqueCombinations.push({
              study: item.study,
              signatureSetName: item.signatureSetName,
              strategy: item.strategy
            });
          }
        });
        
        return uniqueCombinations;
      },
    }),
    satsData: builder.query({
      query: (params) => ({
        url: 'signature_etiology',
        params: { ...params, limit: 1000000 },
      }),
      transformResponse: (data, meta, args) => {
        if (!data || data.length === 0) {
          return { traces: [], layout: {}, config: {} };
        }

        // Transform the etiology data into SATS format
        const transformedData = transformEtiologyDataToSATS(data);
        
        // Generate the SATS plot
        return SATSSignaturePresence(transformedData);
      },
    }),
    satsDataBySignature: builder.query({
      query: (params) => ({
        url: 'mutational_signature',
        params: { 
          signatureName: params.signatureName,
          signatureSetName: params.signatureSetName,
          profile: 'SBS', // Default profile for STS
          matrix: '96'     // Default matrix for STS
        },
      }),
      transformResponse: (data, meta, args) => {
        if (!data || data.length === 0) {
          return { traces: [], layout: {}, config: {} };
        }

        // Transform the mutational signature data into SATS format
        const transformedData = transformEtiologyDataToSATS(data);
        
        // Generate the SATS plot
        return SATSSignaturePresence(transformedData);
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
        SBS: signatureName,
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

export const { useSatsOptionsQuery, useSatsDataQuery, useSatsDataBySignatureQuery, useSatsEtiologyLookupQuery } = satsApiSlice;

// Export transformation functions for reuse
export { transformEtiologyDataToSATS, transformSignatureActivityToSATS };
