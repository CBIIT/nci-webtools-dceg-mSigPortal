import { refittingApiSlice } from '../../../services/store/rootApi';

export const inputFormApiSlice = refittingApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    upload: builder.mutation({
      query: (body) => ({
        url: `uploadRefitting/${crypto.randomUUID()}`,
        method: 'POST',
        body,
      }),
    }),

    submit: builder.mutation({
      query: (body) => ({
        url: `submitRefitting/${body.id}`,
        method: 'POST',
        body,
      }),
    }),

    refresh: builder.query({
      query: (id) => ({
        url: `refreshRefitting/${id}`,
      }),
    }),

    status: builder.query({
      query: (id) => ({
        url: `data/refitting/${id}/status.json`,
      }),
    }),

    params: builder.query({
      query: (id) => ({
        url: `data/refitting/${id}/params.json`,
      }),
    }),

    manifest: builder.query({
      query: (id) => ({
        url: `data/refitting/${id}/manifest.json`,
      }),
    }),

    multiJobStatus: builder.query({
      query: (body) => ({
        url: `refreshRefittingMulti`,
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
              ...status,
            };
          });
      },
    }),

    results: builder.query({
      query: (id) => ({
        url: `data/refitting/${id}/results.json`,
      }),
    }),

    downloadOutput: builder.query({
      query: (id) => ({
        url: `downloadRefittingOutput/${id}`,
      }),
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
  useResultsQuery,
  useDownloadOutputQuery,
} = inputFormApiSlice;
