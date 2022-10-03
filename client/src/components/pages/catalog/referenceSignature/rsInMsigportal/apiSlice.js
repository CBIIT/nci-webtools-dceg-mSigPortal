import { catalogApiSlice } from '../../../../../services/store/rootApi';
import { groupBy } from 'lodash';
import RsInMsigportal from '../../../../controls/plotly/rsInMsigportal/rsInMsigportal';

export const rsInMsigportalApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    rsInMsigportalData: builder.query({
      query: (params) => ({
        url: 'mutational_signature_summary',
        params,
      }),
      transformResponse: (data, meta, args) => {
        console.log(args);
        console.log(data);
        return RsInMsigportal(data);
      },
    }),
  }),
});

export const { useRsInMsigportalDataQuery } = rsInMsigportalApiSlice;
