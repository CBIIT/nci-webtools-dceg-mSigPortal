import { associationApiSlice } from '../../../../services/store/rootApi';

export const publicFormApiSlice = associationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    associationPublic: builder.mutation({
      query: (params) => ({
        url: 'signature_association_options',
        params,
      }),
    }),
  }),
});

export const { useAssociationPublicMutation } = publicFormApiSlice;
