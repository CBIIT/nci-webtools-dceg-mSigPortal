import { explorationApiSlice } from '../../../../services/store/rootApi';

export const userFormApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    exposureUpload: builder.mutation({
      query: (formData) => ({
        url: 'upload',
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
    explorationUser: builder.mutation({
      query: (params) => ({
        url: 'explorationWrapper',
        type: 'POST',
        body: params,
      }),
    }),
  }),
});

export const {
  useExposureUploadMutation,
  useSubmitExplorationMutation,
  useExplorationUserMutation,
} = userFormApiSlice;
