import { apiSlice } from '../../../../services/apiSlice';

export const cosineSimilarityApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    cosineWithin: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    cosineReference: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    cosinePublic: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    cosineSignatureSets: builder.query({
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
  useCosineSignatureSetsQuery,
} = cosineSimilarityApiSlice;
