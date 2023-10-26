import { explorationApiSlice } from '../../../../services/store/rootApi';

export const userFormApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    uploadExploration: builder.mutation({
      query: (formData) => ({
        url: `upload/${crypto.randomUUID()}`,
        method: 'POST',
        body: formData,
      }),
    }),
    submitExploration: builder.mutation({
      query: (data) => ({
        url: `submitExploration/${data.id}`,
        method: 'POST',
        body: data,
      }),
    }),
    explorationWrapper: builder.mutation({
      query: (params) => ({
        url: 'explorationWrapper',
        type: 'POST',
        body: params,
      }),
    }),
  }),
});

export const {
  useUploadExplorationMutation,
  useSubmitExplorationMutation,
  useExplorationWrapperMutation,
} = userFormApiSlice;
