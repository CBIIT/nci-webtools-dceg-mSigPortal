import { extractionApiSlice } from '../../../../services/store/rootApi';

export const inputFormApiSlice = extractionApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    signatureOptions: builder.query({
      query: (params) => ({
        url: 'mutational_signature_options',
        params,
      }),
    }),
  }),
});

export const { useSignatureOptionsQuery } = inputFormApiSlice;
