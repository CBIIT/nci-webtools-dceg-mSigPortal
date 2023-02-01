import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsAssociation from '../../../controls/plotly/msAssociation/msAssociation';
import { groupBy } from 'lodash';

export const msAssociationApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msAssociationOptions: builder.query({
      query: (params) => ({
        url: 'signature_activity',
        params: { ...params },
      }),
      transformResponse: (data) => {
        console.log(data);
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
  }),
});

export const { useMsAssociationOptionsQuery } = msAssociationApiSlice;
