import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsAssociation from '../../../controls/plotly/msAssociation/msAssociation';

export const msAssociationApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msAssociationOptions: builder.query({
      query: (params) => ({
        url: 'signature_activity',
        params: { ...params },
      }),

      transformResponse: (data, meta, arg) => {
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
                id: arg.userId,
              }))
          : [];
      },
    }),
    msAssociation: builder.query({
      query: (params) => ({
        url: 'signature_activity',
        params: params.userId
          ? {
              userId: params.userId,
              signatureName: params.signatureName,
            }
          : {
              study: params.study,
              strategy: params.strategy,
              signatureSetName: params.signatureSetName,
              cancer: params.cancer,
              signatureName: params.signatureName,
            },
      }),
      transformResponse: (data, meta, arg) => {
        const { signatureName, both } = arg;
        const [signatureName1, signatureName2] = signatureName.split(';');
        return MsAssociation(data, signatureName1, signatureName2, both);
      },
    }),
    msAssociation2Source: builder.query({
      async queryFn(args, api, extraOptions, fetchWithBQ) {
        try {
          const [params1, params2] = args;
          const [{ data: data1 }, { data: data2 }] = await Promise.all([
            fetchWithBQ({
              url: '/signature_activity',
              params: {
                userId: params1.userId,
                signatureName: params1.signatureName,
              },
            }),
            fetchWithBQ({
              url: '/signature_activity',
              params: {
                userId: params2.userId,
                signatureName: params2.signatureName,
              },
            }),
          ]);

          const plot = MsAssociation(
            [...data1, ...data2],
            params1.signatureName,
            params2.signatureName,
            params1.both
          );

          return { data: plot };
        } catch (error) {
          console.log(error);
          return { error };
        }
      },
    }),
  }),
});

export const {
  useMsAssociationOptionsQuery,
  useMsAssociationQuery,
  useMsAssociation2SourceQuery,
} = msAssociationApiSlice;
