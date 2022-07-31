import { createApi, fetchBaseQuery } from '@reduxjs/toolkit/query/react';

export const optionsApiSlice = createApi({
  reducerPath: 'optionsApi',
  baseQuery: fetchBaseQuery({ baseUrl: 'web' }),
  endpoints: (builder) => ({
    visualizationOptions: builder.query({
      query: (_) => ({ url: 'seqmatrixOptions' }),
    }),
    explorationOptions: builder.query({
      query: (params) => ({
        url: 'explorationOptions',
        params,
      }),
    }),
  }),
});

export const { useVisualizationOptionsQuery, useExplorationOptionsQuery } =
  optionsApiSlice;

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
  baseQuery: fetchBaseQuery({ baseUrl: '' }),
  endpoints: () => ({}),
});

export const resetOptionsApi = optionsApiSlice.util.resetApiState();
export const resetVisualizationApi = visualizationApiSlice.util.resetApiState();
export const resetExplorationApi = explorationApiSlice.util.resetApiState();
export const resetAssociationApi = associationApiSlice.util.resetApiState();