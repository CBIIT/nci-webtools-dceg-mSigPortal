import { visualizationApiSlice } from '@/services/store/rootApi';

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
      query: ({ strategy, ...params }) => ({
        url: 'mutational_signature_options',
        params,
      }),
      transformResponse: (data, meta, arg) =>
        [
          ...new Set(
            data
              // a de novo set is only valid for the strategy its file was built from
              .filter(
                (e) => e.study === 'Reference' || e.strategy === arg.strategy
              )
              .map((e) => e.signatureSetName)
          ),
        ]
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
