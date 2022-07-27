import { visualizationApiSlice } from '../../../../services/store/rootApi';

export const publicFormApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    publicMatrix: builder.mutation({
      query: (params) => ({
        url: 'seqmatrix',
        params,
      }),
    }),
  }),
});

export const { usePublicMatrixMutation } = publicFormApiSlice;
