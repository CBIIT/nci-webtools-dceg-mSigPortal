import { extractionApiSlice } from '../../../services/store/rootApi';

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
      async queryFn(arg, queryApi, extraOptions, fetchWithBQ) {
        const { id, context_type, signatureMapFile, decomposedSignatureFile } =
          arg;
        const folder = `data/output/${id}/${context_type}/Suggested_Solution/COSMIC_${context_type}_Decomposed_Solution`;
        try {
          const [signatureMap, signatures] = await Promise.all([
            fetchWithBQ(`${folder}/${signatureMapFile}`),
            fetchWithBQ(`${folder}/Signatures/${decomposedSignatureFile}`),
          ]);

          console.log(signatureMap);
          return {
            data: {
              signatureMap: signatureMap.data,
              signatures: signatures.data,
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
  useSignatureMapQuery,
} = inputFormApiSlice;
