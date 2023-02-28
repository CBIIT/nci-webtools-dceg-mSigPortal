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
      transformResponse: (data) => {
        return data.map((job) => {
          const { status, params, manifest } = job;
          if (status)
            return {
              jobName: params.jobName,
              status: status.status,
              id: status.id,
              submittedAt: status.submittedAt,
            };
        });
      },
    }),

    signatureMap: builder.query({
      async queryFn(
        { plotParams, tableParams },
        queryApi,
        extraOptions,
        fetchWithBQ
      ) {
        try {
          const { data: signatureData } = await fetchWithBQ({
            url: `mutational_signature`,
            params: plotParams,
          });
          const { data: tableData } = await fetchWithBQ(
            `data/output/${tableParams.id}/${tableParams.context_type}/Suggested_Solution/COSMIC_${tableParams.context_type}_Decomposed_Solution/${tableParams.signatureMap}`
          );

          const columns = Object.keys(tableData[0]).map((e) => ({
            Header: e,
            accessor: e,
          }));

          // parse signature distribution
          const regex = new RegExp(/Signature\s(\w+)\s\((\d+.\d+)%\)/g);
          const distribution = [
            ...tableData[0]['Global NMF Signatures'].matchAll(regex),
          ].reduce((obj, e) => {
            const [_, key, value] = e;
            return { ...obj, [key]: parseFloat(value) / 100 };
          }, {});
          const decomposedSignatureNames = Object.keys(distribution);
          const denovoSignature = 'SBS96A';

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

          const plots = [denovo, reconstructed, ...decomposed].map((e) =>
            SBS96(e)
          );
          const multi = plots.reduce(
            (obj, e, i) => {
              const count = plots.length - 1;
              const index = i === count ? '' : count - i;
              return {
                traces: [
                  ...obj.traces,
                  ...e.traces.map((trace) => ({
                    ...trace,

                    xaxis: 'x' + index,
                    yaxis: 'y' + index,
                  })),
                ],
                layout: {
                  // ...obj.layout,
                  // ...e.layout,
                  grid: { rows: count, col: 1, pattern: 'independent' },
                },
              };
            },
            {
              traces: [],
              layout: {},
            }
          );

          return {
            data: { plots, table: { data: tableData, columns }, multi },
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
  useSignatureMapQuery,
} = inputFormApiSlice;
