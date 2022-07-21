import { apiSlice } from '../../../../services/apiSlice';

export const publicFormApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    visualizationOptions: builder.query({
      query: (_) => ({ url: 'visualizationOptions' }),
    }),
    visualizationSamples: builder.mutation({
      query: (params) => ({
        url: 'visualizationSamples',
        params,
      }),
    }),
  }),
});

export const { useVisualizationOptionsQuery, useVisualizationSamplesMutation } =
  publicFormApiSlice;
