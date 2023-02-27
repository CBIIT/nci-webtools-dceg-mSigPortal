import { explorationApiSlice } from '../../../../services/store/rootApi';

export const publicFormApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    explorationPublic: builder.mutation({
      query: (params) => ({
        url: 'explorationWrapper',
        type: 'POST',
        body: params,
      }),
    }),
  }),
});

export const { useExplorationPublicMutation } = publicFormApiSlice;
