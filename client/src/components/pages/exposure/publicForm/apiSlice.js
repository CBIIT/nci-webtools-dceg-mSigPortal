import { apiSlice } from '../../../../services/apiSlice';

export const publicFormApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    explorationOptions: builder.query({
      query: (params) => ({
        url: 'explorationOptions',
        params,
      }),
    }),
    explorationSamples: builder.mutation({
      query: (params) => ({
        url: 'explorationSamples',
        params,
      }),
    }),

    explorationPublic: builder.mutation({
      query: (params) => ({
        url: 'explorationWrapper',
        type: 'POST',
        body: params,
      }),
    }),
  }),
});

export const {
  useExplorationOptionsQuery,
  useExplorationSamplesMutation,
  useExplorationPublicMutation,
} = publicFormApiSlice;
