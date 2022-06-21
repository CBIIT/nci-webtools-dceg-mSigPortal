import { createApi, fetchBaseQuery } from '@reduxjs/toolkit/query/react';

export const apiQueries = createApi({
  reducerPath: 'signatureApi',
  baseQuery: fetchBaseQuery({ baseUrl: 'web' }),
  endpoints: (builder) => ({
    signature: builder.query({
      query: (params) => ({ url: 'signature', params }),
    }),
  }),
});

export const { useSignatureQuery } = apiQueries;
