import { extractionApiSlice } from '../../../../services/store/rootApi';

export const inputFormApiSlice = extractionApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    upload: builder.mutation({
      query: (body) => ({
        url: 'upload',
        method: 'POST',
        body,
      }),
    }),
    submit: builder.mutation({
      query: (body) => ({
        url: 'submit',
        method: 'POST',
        body,
      }),
    }),
  }),
});

export const { useUploadMutation, useSubmitMutation } = inputFormApiSlice;
