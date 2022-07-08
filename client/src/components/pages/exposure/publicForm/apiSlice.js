import { apiSlice } from '../../../../services/apiSlice';

export const publicFormApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    getExplorationOptions: builder.query({
      query: (params) => ({
        url: 'explorationOptions',
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
  useGetExplorationOptionsQuery,
  useExplorationPublicMutation,
} = publicFormApiSlice;
