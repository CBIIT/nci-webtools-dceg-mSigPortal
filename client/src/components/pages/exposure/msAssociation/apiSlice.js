import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsAssociation from '../../../controls/plotly/msAssociation/msAssociation';

export const msAssociationApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msAssociationOptions: builder.query({
      query: (params) => ({
        url: 'signature_activity',
        params: { ...params },
      }),

      transformResponse: (data) => {
        return data
          ? [...new Set(data.map((e) => e.signatureName))]
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
    msAssociation: builder.query({
      query: (params) => ({
        url: 'signature_activity',
        params: {
          study: params.study,
          strategy: params.strategy,
          signatureSetName: params.signatureSetName,
          cancer: params.cancer,
          signatureName: params.signatureName,
        },
      }),
      transformResponse: (data, meta, arg) => {
        console.log(arg);

        return MsAssociation(data, arg);
      },
    }),
  }),
});

export const { useMsAssociationOptionsQuery, useMsAssociationQuery } =
  msAssociationApiSlice;
