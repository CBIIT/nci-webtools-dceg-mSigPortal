import { createApi, fetchBaseQuery } from '@reduxjs/toolkit/query/react';

export const apiSlice = createApi({
  reducerPath: 'api',
  baseQuery: fetchBaseQuery({ baseUrl: 'web' }),
  endpoints: (builder) => ({
    seqmatrix: builder.query({
      query: (params) => ({ url: 'seqmatrix', params }),
    }),
    exposure: builder.query({
      query: (params) => ({ url: 'exposure', params }),
    }),
    signature: builder.query({
      query: (params) => ({ url: 'signature', params }),
    }),
  }),
});

export const { useSeqmatrixQuery, useExposureQuery, useSignatureQuery } =
  apiSlice;
