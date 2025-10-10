import { refittingApiSlice } from '../../../services/store/rootApi';

export const inputFormApiSlice = refittingApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    // New SBS refitting endpoint
    submitRefitting: builder.mutation({
      query: (formData) => ({
        url: `refitting/sbs`,
        method: 'POST',
        body: formData,
      }),
    }),

    // Check job status
    refittingStatus: builder.query({
      query: (jobId) => ({
        url: `refitting/status/${jobId}`,
      }),
    }),

    // Download results
    downloadResults: builder.query({
      query: ({ jobId, filename }) => ({
        url: `refitting/download/${jobId}/${filename}`,
      }),
    }),

    // Legacy endpoints (keeping for backward compatibility)
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
  useSubmitRefittingMutation,
  useRefittingStatusQuery,
  useDownloadResultsQuery,
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
