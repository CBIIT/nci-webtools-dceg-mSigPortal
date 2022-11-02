import { catalogApiSlice } from '../../../../services/store/rootApi';
import { groupBy } from 'lodash';
import TMB from '../../../controls/plotly/tmb/tmb';
import SBS96 from '../../../controls/plotly/mutationalProfiles/sbs96';
import SBS192 from '../../../controls/plotly/mutationalProfiles/sbs192';
import DBS78 from '../../../controls/plotly/mutationalProfiles/dbs78';
import ID83 from '../../../controls/plotly/mutationalProfiles/id83';
import CN48 from '../../../controls/plotly/mutationalProfiles/cn48';
import ID29 from '../../../controls/plotly/mutationalProfiles/id29';
import RS32 from '../../../controls/plotly/mutationalProfiles/rs32';

export const etiologyApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    etiologyOptions: builder.query({
      query: (params) => ({
        url: 'signature_etiology_options',
        params,
      }),
    }),
    etiologyDistribtuion: builder.query({
      query: (params) => ({
        url: 'signature_etiology',
        params,
      }),
      transformResponse: (data, _, arg) => {
        if (data.length) {
          const groupByCancer = groupBy(data, 'cancer');
          const transform = Object.entries(groupByCancer)
            .map(([cancer, data]) => {
              const samples = Object.values(groupBy(data, 'sample'))
                .map((e) => ({
                  sample: e[0].sample,
                  burden: e[0].burden,
                  mutations: e[0].mutations,
                }))
                .sort((a, b) => a.burden - b.burden);

              const burdens = samples
                .filter((e) => e.burden)
                .map((e) => e.burden);

              const medianBurden =
                burdens.length % 2 == 0
                  ? (burdens[burdens.length / 2] +
                      burdens[burdens.length / 2 - 1]) /
                    2
                  : burdens[Math.floor(burdens.length / 2)];

              return {
                cancer,
                samples,
                medianBurden,
                sampleSize: data[0].sampleSize,
                signatureSize: data[0].signatureSize,
              };
            })
            .filter((e) => e.medianBurden)
            .sort((a, b) => a.medianBurden - b.medianBurden);

          return TMB(transform, '', arg.signatureName);
        } else {
          return false;
        }
      },
    }),
    etiologySignature: builder.query({
      query: (params) => ({
        url: 'mutational_signature',
        params,
      }),
      transformResponse: (data, _, params) => {
        const profileMatrix = params.profile + params.matrix;
        if (profileMatrix == 'SBS96') return SBS96(data);
        else if (profileMatrix == 'SBS192') return SBS192(data);
        else if (profileMatrix == 'DBS78') return DBS78(data);
        else if (profileMatrix == 'ID83') return ID83(data);
        else if (profileMatrix == 'CN48') return CN48(data);
        else if (profileMatrix == 'ID29') return ID29(data);
        else if (profileMatrix == 'RS32') return RS32(data);
        else return false;
      },
    }),
    etiologyOrganTable: builder.query({
      query: (params) => ({
        url: 'signature_etiology_organ',
        params,
      }),
      transformResponse: (data) => {
        // reducer for creating table columns from objects
        const reducer = (acc, column) => [
          ...acc,
          {
            Header: column,
            id: column,
            accessor: (a) => a[column],
          },
        ];
        return {
          data,
          columns: [...new Set(...data.map((e) => Object.keys(e)))].reduce(
            reducer,
            []
          ),
          cohortOptions: [...new Set(data.map((e) => e.cohort))],
        };
      },
    }),
    thumbnails: builder.query({
      query: (body) => ({
        url: 'getImageS3Batch',
        method: 'POST',
        body,
      }),
    }),
  }),
});

export const {
  useEtiologyOptionsQuery,
  useEtiologyDistribtuionQuery,
  useEtiologySignatureQuery,
  useEtiologyOrganTableQuery,
  useThumbnailsQuery,
} = etiologyApiSlice;
