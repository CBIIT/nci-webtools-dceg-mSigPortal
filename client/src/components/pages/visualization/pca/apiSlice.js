import { visualizationApiSlice } from '@/services/store/rootApi';

export const pcaApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    PcaWithin: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    PcaPublic: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    pcaSignatureSets: builder.query({
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
            a.localeCompare(b, undefined, { sensitivity: 'base' })
          )
          .map((e) => ({ label: e, value: e })),
    }),
  }),
});

export const {
  usePcaWithinQuery,
  usePcaPublicQuery,
  usePcaSignatureSetsQuery,
} = pcaApiSlice;
