import { visualizationApiSlice } from '../../../../services/store/rootApi';

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
      query: (params) => ({ url: 'mutational_signature', params }),
      transformResponse: (data) =>
        [...new Set(data.map((e) => e.signatureSetName))]
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
