import { apiSlice } from '../../../../services/apiSlice';

export const cosineSimilarityApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    CosineWithin: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    CosineReference: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    CosinePublic: builder.query({
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
  }),
});

export const {
  useCosineWithinQuery,
  useCosineReferenceQuery,
  useCosinePublicQuery,
  useSignatureSetsQuery,
} = cosineSimilarityApiSlice;
