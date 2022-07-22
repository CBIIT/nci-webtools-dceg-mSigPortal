import { explorationApiSlice } from '../../../../services/store/rootApi';

export const publicFormApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
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

export const { useExplorationSamplesMutation, useExplorationPublicMutation } =
  publicFormApiSlice;
