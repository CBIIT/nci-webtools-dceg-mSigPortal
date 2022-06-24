import { apiSlice } from '../../../services/apiSlice';

export const vissualizationApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    calculatePublic: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),

    calculateUser: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
  }),
});

export const { useCalculatePublicQuery, useCalculateUserQuery } =
  vissualizationApiSlice;
