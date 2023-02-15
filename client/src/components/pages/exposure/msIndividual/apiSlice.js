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
        //return false;
        return MsIndividual(data, arg);
      },
    }),
  }),
});

export const { useMsIndividualOptionQuery } = msIndividualApiSlice;
