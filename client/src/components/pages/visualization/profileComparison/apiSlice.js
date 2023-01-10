import { visualizationApiSlice } from '../../../../services/store/rootApi';
import sbs96 from '../../../controls/plotly/profileComparision/sbs96';
import sbs192 from '../../../controls/plotly/profileComparision/sbs192';
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
        // get samples from sample array _arg
        const [sample1, sample2] = arg.sample.split(';');
        // filter api data by samples
        const userData = data.filter(
          (e) => e.sample == sample1 || e.signatureName == sample1
        );
        const publicData = data.filter(
          (e) => e.sample == sample2 || e.signatureName == sample2
        );

        if (arg.profile === 'SBS') {
          return sbs96(userData, publicData);
        } else if (arg.profile === 'DBS') {
          return dbs78(userData, publicData);
        } else if (arg.profile == 'ID') {
          return id83(userData, publicData);
        } else throw Error(`Profile ${arg.profile} is not supported`);
      },
    }),

    profileComparisonReference: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        try {
          console.log(_arg.params_signature);
          console.log(_arg.params_signature_scalar);

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
          const scalarArray = Object.values(_arg.params_signature_scalar)[0];
          const signatureArray = Object.values(_arg.params_signature_scalar)[1];

          console.log(res);
          console.log(signatureData);
          console.log(_arg.params_signature_scalar);
          console.log(scalarArray);
          console.log(signatureArray);

          let result = [];

          // for (var i = 0; i < )

          if (_arg.params_spectrum.profile === 'SBS') {
            return {
              data: sbs96(spectrumData, signatureData),
            };
          } else if (_arg.params_spectrum.profile === 'DBS') {
            return {
              data: dbs78(spectrumData, signatureData),
            };
          } else {
            return {
              data: id83(spectrumData, signatureData),
            };
          }
        } catch (error) {
          return { error };
        }
      },
    }),

    profileComparisonPublic: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        const { profile, matrix } = _arg.userParams;

        try {
          const [{ data: userData }, { data: publicData }] = await Promise.all([
            fetchWithBQ(
              '/mutational_spectrum?' + new URLSearchParams(_arg.userParams)
            ),
            fetchWithBQ(
              '/mutational_spectrum?' + new URLSearchParams(_arg.publicParams)
            ),
          ]);

          if (profile == 'SBS' && matrix == '96')
            return { data: sbs96(userData, publicData) };
          else if (profile == 'SBS' && matrix == '192')
            return { data: sbs192(userData, publicData) };
          else if (profile == 'DBS')
            return { data: dbs78(userData, publicData) };
          else if (profile == 'ID') return { data: id83(userData, publicData) };
          else
            throw Error(
              `Profile/Matrix: ${profile}/${matrix} is not supported`
            );
        } catch (error) {
          return { error: error };
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
