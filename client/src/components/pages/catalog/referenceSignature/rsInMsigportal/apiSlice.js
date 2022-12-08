import { catalogApiSlice } from '../../../../../services/store/rootApi';
import RsInMsigportal from '../../../../controls/plotly/rsInMsigportal/rsInMsigportal';

export const rsInMsigportalApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    rsInMsigportalData: builder.query({
      query: (params) => ({
        url: 'mutational_signature_summary',
        params,
      }),
      transformResponse: (data, meta, args) => {
        return RsInMsigportal(data);
      },
    }),
  }),
});

export const { useRsInMsigportalDataQuery } = rsInMsigportalApiSlice;
