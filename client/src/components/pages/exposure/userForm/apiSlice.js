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
    exposureUser: builder.mutation({
      query: (data) => ({
        url: 'getSignaturesUser',
        method: 'POST',
        body: data,
      }),
    }),
  }),
});

export const { useExposureUploadMutation, useExposureUserMutation } =
  userFormApiSlice;
