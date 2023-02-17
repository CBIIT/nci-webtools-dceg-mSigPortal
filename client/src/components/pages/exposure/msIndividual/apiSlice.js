import { explorationApiSlice } from '../../../../services/store/rootApi';
import msIndividual_sbs96 from '../../../controls/plotly/msIndividual/sbs96';

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
    msIndividual: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        try {
          const res = await Promise.all([
            fetchWithBQ(
              '/signature_activity?' + new URLSearchParams(_arg.params_activity)
            ), //exposure
            fetchWithBQ(
              '/mutational_signature?' +
                new URLSearchParams(_arg.params_signature)
            ), //signature
            fetchWithBQ(
              '/mutational_spectrum?' +
                new URLSearchParams(_arg.params_spectrum)
            ), //seqmatrix
          ]);

          console.log(res);
          console.log(_arg);
          //return MsLandscape(res, _arg);
          return { data: msIndividual_sbs96(res, _arg) };
        } catch (error) {
          return { error };
        }
      },
    }),
  }),
});

export const { useMsIndividualOptionQuery, useMsIndividualQuery } =
  msIndividualApiSlice;
