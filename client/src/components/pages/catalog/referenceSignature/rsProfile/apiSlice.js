import { catalogApiSlice } from '../../../../services/store/rootApi';
import { groupBy } from 'lodash';

export const rsProfileApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    rsProfile: builder.query({
      query: (params) => ({
        url: 'mutational_signature',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        //   if (arg.profile === 'SBS') {
        //     return SBS96(samples, sample1, sample2, 'samples');
        //   } else if (arg.profile === 'DBS') {
        //     return DBS78(samples, sample1, sample2, 'samples');
        //   } else {
        //     return ID83(samples, sample1, sample2, 'samples');
        //   }
        console.log(data);
      },
    }),
  }),
});

export const { useRsProfileQuery } = rsProfileApiSlice;
