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
  }),
});

// Transform function to convert etiology API data to SATS format
function transformEtiologyDataToSATS(etiologyData) {
  if (!etiologyData || etiologyData.length === 0) {
    return { tmbData: [], dotData: [] };
  }

  // Group data by cancer type first
  const groupedByCancer = groupBy(etiologyData, 'cancer');
  
  const satsData = [];
  let cancerTypeNum = 1;
  
  Object.entries(groupedByCancer).forEach(([cancerType, cancerData]) => {
    // Group by signature for this cancer type
    const groupedBySignature = groupBy(cancerData, 'signatureName');
    
    // Calculate total mutations across all samples for this cancer type
    const totalMutations = cancerData.reduce((sum, item) => sum + (item.mutations || 0), 0);
    const totalSamples = [...new Set(cancerData.map(item => item.sample))].length; // unique samples
    
    // Calculate TMB_all (total mutational burden for this cancer type)
    const tmbAll = totalSamples > 0 ? totalMutations / totalSamples : 0;
    
    Object.entries(groupedBySignature).forEach(([signatureName, signatureData]) => {
      // Calculate metrics for this cancer-signature combination
      const totalExposure = signatureData.reduce((sum, item) => sum + (item.exposure || 0), 0);
      
      // Count samples with detectable signature (exposure > 0)
      const samplesWithSignature = signatureData.filter(item => (item.exposure || 0) > 0).length;
      
      // Calculate presence (proportion of samples with this signature)
      const presence = totalSamples > 0 ? samplesWithSignature / totalSamples : 0;
      
      // Calculate proportion (signature's contribution to total mutations)
      const proportion = totalMutations > 0 ? totalExposure / totalMutations : 0;
      
      satsData.push({
        CancerType: cancerType,
        SBS: signatureName,
        Count: totalExposure,
        Proportion: proportion,
        N: totalSamples,
        Presence: presence,
        TMB_all: tmbAll,
        Label: totalSamples,
        CancerType_num: cancerTypeNum
      });
    });
    
    cancerTypeNum++;
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

export const { useSatsOptionsQuery, useSatsDataQuery } = satsApiSlice;

// Export transformation functions for reuse
export { transformEtiologyDataToSATS, transformSignatureActivityToSATS };
