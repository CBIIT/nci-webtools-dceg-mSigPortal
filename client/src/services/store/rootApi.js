import { createApi, fetchBaseQuery } from '@reduxjs/toolkit/query/react';

export const optionsApiSlice = createApi({
  reducerPath: 'optionsApi',
  baseQuery: fetchBaseQuery({ baseUrl: 'web' }),
  endpoints: (builder) => ({
    seqmatrixOptions: builder.query({
      query: (params) => ({ url: 'seqmatrixOptions', params }),
    }),
    exposureOptions: builder.query({
      query: (params) => ({ url: 'exposureOptions', params }),
    }),
    signatureOptions: builder.query({
      query: (params) => ({ url: 'signatureOptions', params }),
    }),
    associationOptions: builder.query({
      query: (params) => ({ url: 'associationOptions', params }),
    }),
  }),
});

export const {
  useSeqmatrixOptionsQuery,
  useExposureOptionsQuery,
  useSignatureOptionsQuery,
  useAssociationOptionsQuery,
} = optionsApiSlice;

export const visualizationApiSlice = createApi({
  reducerPath: 'visualizationApi',
  baseQuery: fetchBaseQuery({ baseUrl: 'web' }),
  endpoints: () => ({}),
});

export const explorationApiSlice = createApi({
  reducerPath: 'explorationApi',
  baseQuery: fetchBaseQuery({ baseUrl: 'web' }),
  endpoints: () => ({}),
});

export const associationApiSlice = createApi({
  reducerPath: 'associationApi',
  baseQuery: fetchBaseQuery({ baseUrl: 'web' }),
  endpoints: () => ({}),
});

export const catalogApiSlice = createApi({
  reducerPath: 'catalogApi',
  baseQuery: fetchBaseQuery({ baseUrl: 'web' }),
  endpoints: () => ({}),
});

export const resetOptionsApi = optionsApiSlice.util.resetApiState();
export const resetVisualizationApi = visualizationApiSlice.util.resetApiState();
export const resetExplorationApi = explorationApiSlice.util.resetApiState();
export const resetAssociationApi = associationApiSlice.util.resetApiState();
