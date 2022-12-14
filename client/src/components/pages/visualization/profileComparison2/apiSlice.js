import { visualizationApiSlice } from '../../../../services/store/rootApi';
import sbs96 from '../../../controls/plotly/profileComparision/sbs96';
import dbs78 from '../../../controls/plotly/profileComparision/dbs78';
import id83 from '../../../controls/plotly/profileComparision/id83';

export const profilerSummaryApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    profileComparisonWithin: builder.query({
      query: (params) => ({
        url: 'mutational_spectrum',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        const samples = arg.sample.split(';');

        if (arg.profile === 'SBS') {
          return sbs96(samples, data);
        } else if (arg.profile === 'DBS') {
          return dbs78(samples, data);
        } else {
          return id83(samples, data);
        }
      },
    }),

    profileComparisonReference: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        try {
          const res = await Promise.all([
            fetchWithBQ(
              '/mutational_spectrum?' +
                new URLSearchParams(_arg.params_spectrum)
            ),
            fetchWithBQ(
              '/mutational_signature?' +
                new URLSearchParams(_arg.params_signature)
            ),
          ]);

          const spectrumData = res[0]['data'];
          const signatureData = res[1]['data'];
          const samples = [
            _arg.params_spectrum.sample,
            _arg.params_signature.signatureName,
          ];

          if (_arg.params_spectrum.profile === 'SBS') {
            return {
              data: sbs96(samples, spectrumData, signatureData, 'reference'),
            };
          } else if (_arg.params_spectrum.profile === 'DBS') {
            return {
              data: dbs78(samples, spectrumData, signatureData, 'reference'),
            };
          } else {
            return {
              data: id83(samples, spectrumData, signatureData, 'reference'),
            };
          }
        } catch (error) {
          return { error };
        }
      },
    }),

    profileComparisonPublic: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        try {
          const [userData, publicData] = await Promise.all([
            fetchWithBQ(
              '/mutational_spectrum?' + new URLSearchParams(_arg.userParams)
            ), //seqmatrix
            fetchWithBQ(
              '/mutational_spectrum?' + new URLSearchParams(_arg.publicParams)
            ),
          ]);

          return { data: { userData, publicData } };
        } catch (error) {
          return { error };
        }
      },
    }),

    pcSignatureSets: builder.query({
      query: (params) => ({ url: 'mutational_signature', params }),
      transformResponse: (data) =>
        [...new Set(data.map((e) => e.signatureSetName))]
          .sort((a, b) =>
            a.localeCompare(b, undefined, { sensitivity: 'base' })
          )
          .map((e) => ({ label: e, value: e })),
    }),
    pcSignatureNames: builder.query({
      query: (params) => ({ url: 'mutational_signature', params }),
      transformResponse: (data, meta, arg) =>
        [...new Set(data.map((e) => e.signatureName))].sort((a, b) =>
          a.localeCompare(b, undefined, { sensitivity: 'base' })
        ),
    }),
  }),
});

export const {
  useProfileComparisonWithinQuery,

  useProfileComparisonReferenceQuery,
  useProfileComparisonPublicQuery,
  usePcSignatureSetsQuery,
  usePcSignatureNamesQuery,
} = profilerSummaryApiSlice;
