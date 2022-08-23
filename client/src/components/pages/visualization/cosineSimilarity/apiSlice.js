import { visualizationApiSlice } from '../../../../services/store/rootApi';

export const cosineSimilarityApiSlice = visualizationApiSlice.injectEndpoints({
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
      query: (params) => ({ url: 'mutational_signature', params }),
      transformResponse: (data) =>
        [...new Set(data.map((e) => e.signatureSetName))]
          .sort((a, b) =>
            a.localeCompare(b, undefined, {
              numeric: true,
              sensitivity: 'base',
            })
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
