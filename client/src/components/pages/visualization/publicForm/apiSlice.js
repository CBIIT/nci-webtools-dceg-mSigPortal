import { visualizationApiSlice } from '../../../../services/store/rootApi';

export const publicFormApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    visualizationSamples: builder.mutation({
      query: (params) => ({
        url: 'seqmatrix',
        params,
      }),
    }),
  }),
});

export const { useVisualizationSamplesMutation } = publicFormApiSlice;
