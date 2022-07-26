import { createApi, fetchBaseQuery } from '@reduxjs/toolkit/query/react';

export const optionsApiSlice = createApi({
  reducerPath: 'optionsApi',
  baseQuery: fetchBaseQuery({ baseUrl: '' }),
  endpoints: (builder) => ({
    visualizationOptions: builder.query({
      query: (_) => ({ url: 'visualizationOptions' }),
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
  baseQuery: fetchBaseQuery({ baseUrl: '' }),
  endpoints: () => ({}),
});

export const explorationApiSlice = createApi({
  reducerPath: 'explorationApi',
  baseQuery: fetchBaseQuery({ baseUrl: '' }),
  endpoints: () => ({}),
});

export const resetOptionsApi = optionsApiSlice.util.resetApiState();
export const resetVisualizationApi = visualizationApiSlice.util.resetApiState();
export const resetExplorationApi = explorationApiSlice.util.resetApiState();
