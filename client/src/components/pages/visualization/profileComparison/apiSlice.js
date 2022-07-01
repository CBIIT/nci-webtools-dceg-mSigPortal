import { apiSlice } from '../../../../services/apiSlice';

export const profilerSummaryApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    profileComparisonWithin: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
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
    signatureSets: builder.query({
      query: (params) => ({ url: 'signature', params }),
      transformResponse: (data) =>
        data
          .map((e) => e.signatureSetName)
          .sort((a, b) =>
            a.localeCompare(b, undefined, { sensitivity: 'base' })
          )
          .map((e) => ({ label: e, value: e })),
    }),
    signatureNames: builder.query({
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
  useSignatureSetsQuery,
  useSignatureNamesQuery,
} = profilerSummaryApiSlice;
