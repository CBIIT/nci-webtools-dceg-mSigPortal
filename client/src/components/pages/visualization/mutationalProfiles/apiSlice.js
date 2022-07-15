import { apiSlice } from '../../../../services/apiSlice';

export const mutationalProfilesApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    mutationalProfiles: builder.query({
      query: (params) => ({ url: 'mutationalProfiles', params }),
    }),
  }),
});

export const { useMutationalProfilesQuery } = mutationalProfilesApiSlice;
