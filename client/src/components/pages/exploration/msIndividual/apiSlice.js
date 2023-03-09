import { explorationApiSlice } from '../../../../services/store/rootApi';
import sbs96 from '../../../controls/plotly/profileComparison/sbs96';
import dbs78 from '../../../controls/plotly/profileComparison/dbs78';
import id83 from '../../../controls/plotly/profileComparison/id83';
import msIndividual_rs32 from '../../../controls/plotly/msIndividual/msIndividual_rs32';
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
        const res = await Promise.all([
          fetchWithBQ(
            '/signature_activity?' + new URLSearchParams(_arg.params_activity)
          ), //exposure
          fetchWithBQ(
            '/mutational_signature?' +
              new URLSearchParams(_arg.params_signature)
          ), //signature
          fetchWithBQ(
            '/mutational_spectrum?' + new URLSearchParams(_arg.params_spectrum)
          ), //seqmatrix
        ]);
        let profile;

        if (_arg.params_activity.userId) {
          if (res[0].data.length > 0) {
            if (res[0].data[0].signatureName === res[1].data[0].signatureName) {
              profile = res[0].data[0].signatureName;
            } else {
              profile = 'Data mismatch';
            }
          } else {
            profile = 'No data';
          }
        } else {
          if (res[0].data.length > 0) {
            profile = extractLastWord(_arg.params_activity.signatureSetName);
          } else {
            profile = 'No data';
          }
        }

        if (profile === 'SBS96' || profile.includes('SBS')) {
          return { data: sbs96(res, _arg, 'msIndividual') };
        } else if (profile === 'DBS78' || profile.includes('DBS')) {
          return { data: dbs78(res, _arg, 'msIndividual') };
        } else if (profile === 'ID83' || profile.includes('ID')) {
          return { data: id83(res, _arg, 'msIndividual') };
        } else {
          //return { data: msIndividual_rs32(res, _arg, 'msIndividual') };
          return {
            error: _arg.params_activity.signatureSetName
              ? 'Signature SetName: ' +
                _arg.params_activity.signatureSetName +
                ' is not supported in MS Individual'
              : profile === 'Data mismatch'
              ? 'Data mismatch, please try again with different files'
              : 'No data found, please try again',
          };
        }
      },
    }),
  }),
});

export const { useMsIndividualOptionQuery, useMsIndividualQuery } =
  msIndividualApiSlice;
