import { associationApiSlice } from '../../../../services/store/rootApi';

export const publicFormApiSlice = associationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    associationPublic: builder.mutation({
      query: (params) => ({
        url: 'associationOptions',
        params,
      }),
    }),
  }),
});

export const { useAssociationPublicMutation } = publicFormApiSlice;
