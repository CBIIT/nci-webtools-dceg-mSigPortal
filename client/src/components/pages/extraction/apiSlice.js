import { extractionApiSlice } from '../../../services/store/rootApi';
import { groupBy } from 'lodash';
import SBS96 from '../../controls/plotly/mutationalProfiles/sbs96';
import DBS78 from '../../controls/plotly/mutationalProfiles/dbs78';
import ID83 from '../../controls/plotly/mutationalProfiles/id83';

export const inputFormApiSlice = extractionApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    upload: builder.mutation({
      query: (body) => ({
        url: `upload/${crypto.randomUUID()}`,
        method: 'POST',
        body,
      }),
    }),

    submit: builder.mutation({
      query: (body) => ({
        url: `submitExtraction/${body.id}`,
        method: 'POST',
        body,
      }),
    }),

    refresh: builder.query({
      query: (id) => ({
        url: `refreshExtraction/${id}`,
      }),
    }),

    status: builder.query({
      query: (id) => ({
        url: `data/output/${id}/status.json`,
      }),
    }),

    params: builder.query({
      query: (id) => ({
        url: `data/input/${id}/params.json`,
      }),
    }),

    manifest: builder.query({
      query: (id) => ({
        url: `data/output/${id}/manifest.json`,
      }),
    }),

    multiJobStatus: builder.query({
      query: (body) => ({
        url: `refreshExtractionMulti`,
        method: 'POST',
        body,
      }),
      transformResponse: (data, meta, args) => {
        return data
          .filter((e) => e.status && e.params)
          .map((job) => {
            const { status, params, manifest } = job;
            return {
              jobName: params.jobName,
              status: status.status,
              id: status.id,
              submittedAt: status.submittedAt,
            };
          });
      },
    }),

    signatureMapTable: builder.query({
      query: ({ id, context_type, signatureMap }) => ({
        url: `data/output/${id}/${context_type}/Suggested_Solution/COSMIC_${context_type}_Decomposed_Solution/${signatureMap}`,
      }),
      transformResponse: (data) => {
        const columns = Object.keys(data[0]).map((e) => ({
          Header: e,
          accessor: e,
        }));

        return { data, columns };
      },
    }),

    signatureMapPlots: builder.query({
      async queryFn(params, queryApi, extraOptions, fetchWithBQ) {
        try {
          const profileMatrixMap = { SBS: 96, DBS: 78, ID: 83 };
          const { userId, decompSigString, denovoSigString } = params;

          // parse signature distribution names and proportion
          const distributionRegex = new RegExp(
            /Signature\s(\w+)\s\((\d+.\d+)%\)/g
          );
          const distribution = [
            ...decompSigString.matchAll(distributionRegex),
          ].reduce((obj, e) => {
            const [_, key, value] = e;
            return { ...obj, [key]: parseFloat(value) / 100 };
          }, {});
          const decomposedSignatureNames = Object.keys(distribution);

          // parse denovo signature name
          const profileRegex = /([a-zA-Z]+)/;
          const profile = decomposedSignatureNames[0].match(profileRegex)[1];
          const profileMatrix = profile + profileMatrixMap[profile];
          const denovoSignature = profileMatrix + denovoSigString.slice(-1);

          // query signatures
          const { data: signatureData } = await fetchWithBQ({
            url: `mutational_signature`,
            params: { userId },
          });

          const allSignatures = groupBy(signatureData, (e) => e.signatureName);
          const reconstructed = Object.values(
            groupBy(signatureData, (e) => e.mutationType)
          ).map((data) =>
            data
              .filter((e) => decomposedSignatureNames.includes(e.signatureName))
              .reduce(
                (obj, e) => ({
                  ...obj,
                  ...e,
                  mutations:
                    (obj?.mutations || 0) +
                    e.mutations * distribution[e.signatureName],
                  signatureName: `${denovoSignature} (Reconstructed)`,
                }),
                {}
              )
          );

          const denovo = allSignatures[denovoSignature];

          const decomposed = decomposedSignatureNames.map((e) =>
            allSignatures[e].reduce(
              (array, e) => [
                ...array,
                {
                  ...e,
                  signatureName: `${e.signatureName} (${(
                    distribution[e.signatureName] * 100
                  ).toFixed(2)}%)`,
                },
              ],
              []
            )
          );

          const plots = [
            { title: 'Denovo Signature', data: denovo },
            { title: 'Reconstructed Signature', data: reconstructed },
            ...decomposed.map((e) => ({
              title: 'Decomposed Signature',
              data: e,
            })),
          ].map((e) => {
            switch (profileMatrix) {
              case 'SBS96':
                return SBS96(e.data, e.title);
              case 'DBS78':
                return DBS78(e.data, e.title);
              case 'ID83':
                return ID83(e.data, e.title);
              default:
                throw new Error(`${profileMatrix} is not supported`);
            }
          });

          return {
            data: { plots },
          };
        } catch (error) {
          return { error };
        }
      },
    }),
  }),
});

export const {
  useUploadMutation,
  useSubmitMutation,
  useRefreshQuery,
  useStatusQuery,
  useParamsQuery,
  useManifestQuery,
  useMultiJobStatusQuery,
  useSignatureMapTableQuery,
  useSignatureMapPlotsQuery,
} = inputFormApiSlice;
