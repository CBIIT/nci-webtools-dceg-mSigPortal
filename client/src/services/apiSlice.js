import { createApi, fetchBaseQuery } from '@reduxjs/toolkit/query/react';

export const apiSlice = createApi({
  reducerPath: 'api',
  baseQuery: fetchBaseQuery({ baseUrl: 'web' }),
  endpoints: (builder) => ({
    signature: builder.query({
      query: (params) => ({ url: 'signature', params }),
    }),
    exposure: builder.query({
      query: (params) => ({ url: 'exposure', params }),
    }),
    signatureRefSets: builder.query({
      query: (params) => ({ url: 'signatureRefSets', params }),
    }),
  }),
});

export const { useSignatureQuery, useExposureQuery, useSignatureRefSetsQuery } =
  apiSlice;
