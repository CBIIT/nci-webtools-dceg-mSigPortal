import { explorationApiSlice } from '../../../../services/store/rootApi';
import sbs96 from '../../../controls/plotly/profileComparison/sbs96';
import dbs78 from '../../../controls/plotly/profileComparison/dbs78';
import id83 from '../../../controls/plotly/profileComparison/id83';
import msIndividual_rs32 from '../../../controls/plotly/profileComparison/rs32';
import { extractLastWord } from '../../../controls/utils/utils';

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
          const profile = extractLastWord(
            _arg.params_activity.signatureSetName
          );
          console.log(profile);
          if (profile === 'SBS96') {
            return { data: sbs96(res, _arg, 'msIndividual') };
          } else if (profile === 'DBS78') {
            return { data: dbs78(res, _arg, 'msIndividual') };
          } else if (profile === 'ID83') {
            return { data: id83(res, _arg, 'msIndividual') };
          } else if (profile === 'RS32') {
            return { data: msIndividual_rs32(res, _arg, 'msIndividual') };
          } else
            throw Error(
              `Signature Set Name: ${_arg.params_activity.signatureSetName} is not supported`
            );
        } catch (error) {
          return { error };
        }
      },
    }),
  }),
});

export const { useMsIndividualOptionQuery, useMsIndividualQuery } =
  msIndividualApiSlice;
