import { apiSlice } from '../../../../services/apiSlice';
import MsBurden from '../../../controls/plotly//msBurden/msBurden';

export const msBurdenApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msBurden: builder.query({
      query: (params) => ({
        url: 'tmb',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        return MsBurden(data);
      },
    }),
  }),
});

export const { useMsBurdenQuery } = msBurdenApiSlice;
