import { visualizationApiSlice } from '../../../../services/store/rootApi';
import { parseCSV } from '../../../../services/utils';

export const userFormApiSlice = visualizationApiSlice.injectEndpoints({
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
        url: `submitVisualization/${body.id}`,
        method: 'POST',
        body,
      }),
    }),

    refresh: builder.query({
      async queryFn(id, api, extraOptions, fetchWithBQ) {
        try {
          const [{ data: status }, { data: params }, { data: manifest }] =
            await Promise.all([
              fetchWithBQ(`data/output/${id}/status.json`),
              fetchWithBQ(`data/input/${id}/params.json`),
              fetchWithBQ(`data/output/${id}/manifest.json`),
            ]);
          return { data: { status, params, manifest } };
        } catch (error) {
          return { error };
        }
      },
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

    matrixList: builder.query({
      query: (id) => ({
        url: `data/output/${id}/profilerExtraction/matrix_files_list.txt`,
        responseHandler: (response) => response.text(),
      }),
      transformResponse: async (data) => {
        const parsed = await parseCSV(data);
        const matrixList = parsed.map(
          ({ Profile_Type, Matrix_Size, ...e }) => ({
            ...e,
            profile: Profile_Type,
            matrix: Matrix_Size,
          })
        );
        return matrixList;
      },
    }),

    multiJobStatus: builder.query({
      query: (body) => ({
        url: `refreshVisualizationMulti`,
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

    exampleHeader: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        try {
          const data = await (
            await fetch(
              `https://raw.githubusercontent.com/xtmgah/mSigPortal/master/Input_Format/Header_${_arg}.txt`
            )
          ).text();
          return { data };
        } catch (error) {
          console.log(error);
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
  useMatrixListQuery,
  useMultiJobStatusQuery,
  useExampleHeaderQuery,
} = userFormApiSlice;
