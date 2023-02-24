import { visualizationApiSlice } from '../../../../services/store/rootApi';
import sbs96 from '../../../controls/plotly/profileComparison/sbs96';
import sbs192 from '../../../controls/plotly/profileComparison/sbs192';
import dbs78 from '../../../controls/plotly/profileComparison/dbs78';
import id83 from '../../../controls/plotly/profileComparison/id83';
import { groupBy } from 'lodash';

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
          return sbs96(userData, publicData, 'pc');
        } else if (arg.profile === 'DBS') {
          return dbs78(userData, publicData, 'pc');
        } else if (arg.profile == 'ID') {
          return id83(userData, publicData, 'pc');
        } else throw Error(`Profile ${arg.profile} is not supported`);
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

          // normalize signatureData by taking the average of contribution of selected signatureNames
          const groupByMutation = Object.values(
            groupBy(signatureData, (e) => e.mutationType)
          );
          const normalizeSigData = groupByMutation.map((signatures) => ({
            ...signatures[0],
            signatureName: _arg.params_signature.signatureName,
            scalarSignature: _arg.params_signature_scalar.paramsSignatureScalar,
            contribution:
              signatures.reduce((total, e) => total + e.contribution, 0) /
              signatures.length,
          }));

          if (_arg.params_spectrum.profile === 'SBS') {
            return {
              data: sbs96(spectrumData, normalizeSigData),
            };
          } else if (_arg.params_spectrum.profile === 'DBS') {
            return {
              data: dbs78(spectrumData, normalizeSigData),
            };
          } else {
            return {
              data: id83(spectrumData, normalizeSigData),
            };
          }
        } catch (error) {
          console.error(error);
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
  }),
});

export const {
  useProfileComparisonWithinQuery,
  useProfileComparisonReferenceQuery,
  useProfileComparisonPublicQuery,
} = profilerSummaryApiSlice;
