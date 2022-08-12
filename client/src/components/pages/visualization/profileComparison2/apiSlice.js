import { visualizationApiSlice } from '../../../../services/store/rootApi';
import pcBetweenSamples from '../../../controls/plotly/profileComparision/pcBetweenSamples';
export const profilerSummaryApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    // profileComparisonWithin: builder.query({
    //   query: (params) => ({
    //     url: 'visualizationWrapper',
    //     method: 'POST',
    //     body: params,
    //   }),
    // }),
    profileComparisonWithin: builder.query({
      query: (params) => ({
        url: 'seqmatrix',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        console.log(data);
        return pcBetweenSamples(data, arg);
      },
    }),
    profileComparisonReference: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    profileComparisonPublic: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    pcSignatureSets: builder.query({
      query: (params) => ({ url: 'signature', params }),
      transformResponse: (data) =>
        [...new Set(data.map((e) => e.signatureSetName))]
          .sort((a, b) =>
            a.localeCompare(b, undefined, { sensitivity: 'base' })
          )
          .map((e) => ({ label: e, value: e })),
    }),
    pcSignatureNames: builder.query({
      query: (params) => ({ url: 'signature', params }),
      transformResponse: (data) =>
        [...new Set(data.map((e) => e.signatureName))].sort((a, b) =>
          a.localeCompare(b, undefined, { sensitivity: 'base' })
        ),
    }),
  }),
});

export const {
  useProfileComparisonWithinQuery,
  useProfileComparisonReferenceQuery,
  useProfileComparisonPublicQuery,
  usePcSignatureSetsQuery,
  usePcSignatureNamesQuery,
} = profilerSummaryApiSlice;
