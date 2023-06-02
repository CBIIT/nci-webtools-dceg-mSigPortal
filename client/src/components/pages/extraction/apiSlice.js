import { extractionApiSlice } from '../../../services/store/rootApi';
import { groupBy } from 'lodash';

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

    example: builder.query({
      query: (id) => ({
        url: `extractionExample/${id}`,
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

          let denovo;
          let reconstructed;
          let decomposed;
          let profileMatrix;
          let refSigPlots;
          let denovoPlots;

          async function createPlot(profileMatrix, data, title) {
            try {
              const plotFn = await import(
                `../../controls/plotly/mutationalProfiles/${profileMatrix.toLowerCase()}.js`
              );
              return plotFn.default(data, title);
            } catch (error) {
              console.log(error);
              throw new Error(`${profileMatrix} is not supported`);
            }
          }

          // query signatures
          const { data: signatureData } = await fetchWithBQ({
            url: `mutational_signature`,
            params: { userId },
          });
          const allSignatures = groupBy(signatureData, (e) => e.signatureName);
          if (decompSigString === denovoSigString) {
            const filteredSignatures = {};

            for (const key in allSignatures) {
              if (key.endsWith(denovoSigString.split('-')[1])) {
                filteredSignatures[key] = allSignatures[key];
              }
            }

            profileMatrix = Object.keys(filteredSignatures)[0].substring(
              0,
              Object.keys(filteredSignatures)[0].length - 1
            );

            const groupByMutationType = groupBy(
              Object.values(filteredSignatures)[0],
              'mutationType'
            );

            const firstValues = Object.keys(groupByMutationType).map(
              (key) => groupByMutationType[key][0]
            );

            denovo = firstValues;
            decomposed = firstValues;

            refSigPlots = [];
            denovoPlots = await Promise.all(
              [{ title: 'Denovo Signature', data: denovo }].map(
                async (e) => await createPlot(profileMatrix, e.data, e.title)
              )
            );
          } else {
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
            profileMatrix = profile + profileMatrixMap[profile];
            const denovoSignature = profileMatrix + denovoSigString.slice(-1);

            reconstructed = Object.values(
              groupBy(signatureData, (e) => e.mutationType)
            ).map((data) =>
              data
                .filter((e) =>
                  decomposedSignatureNames.includes(e.signatureName)
                )
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

            denovo = allSignatures[denovoSignature].map((e) => ({
              ...e,
              signatureName: `${denovoSignature} (Original)`,
            }));

            decomposed = decomposedSignatureNames.map((e) =>
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

            refSigPlots = {};
            decomposed.forEach(async (e) => {
              refSigPlots[e[0].signatureName] = await createPlot(
                profileMatrix,
                e,
                'Reference Signature'
              );
            });

            denovoPlots = await Promise.all(
              [
                { title: 'Denovo Signature', data: denovo },
                { title: 'Reconstructed Signature', data: reconstructed },
              ].map((e) => createPlot(profileMatrix, e.data, e.title))
            );
          }

          // parse signature distribution names and proportion
          return {
            data: {
              denovoPlots: denovoPlots,
              refSigPlots: refSigPlots,
            },
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
  useExampleQuery,
} = inputFormApiSlice;
