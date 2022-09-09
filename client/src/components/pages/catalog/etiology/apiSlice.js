import { catalogApiSlice } from '../../../../services/store/rootApi';

export const etiologyApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    etiologyOptions: builder.query({
      query: (params) => ({
        url: 'etiologyOptions',
        params,
      }),
    }),
    thumbnails: builder.query({
      query: (body) => ({
        url: 'getImageS3Batch',
        method: 'POST',
        body,
      }),
    }),
  }),
});

export const { useEtiologyOptionsQuery, useThumbnailsQuery } = etiologyApiSlice;
