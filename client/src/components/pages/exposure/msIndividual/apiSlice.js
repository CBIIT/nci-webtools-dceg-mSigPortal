import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsIndividual from '../../../controls/plotly/msIndividual/msIndividual';

export const msIndividualApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msIndividualOption: builder.query({
      query: (params) => ({
        url: 'signature_activity',
        params: {
          ...params,
        },
      }),
      transformResponse: (data, meta, arg) => {
        console.log('---ms indiidual');
        console.log(arg);
        console.log(data);
        return data
          ? [...new Set(data.map((e) => e.sample))]
              .sort((a, b) =>
                a.localeCompare(b, undefined, {
                  numeric: true,
                  sensitivity: 'base',
                })
              )
              .map((e) => ({
                label: e,
                value: e,
              }))
          : [];
      },
    }),
  }),
});

export const { useMsIndividualOptionQuery } = msIndividualApiSlice;
