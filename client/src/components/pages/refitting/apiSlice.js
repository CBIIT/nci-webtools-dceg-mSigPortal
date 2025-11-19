import { refittingApiSlice } from '../../../services/store/rootApi';

export const inputFormApiSlice = refittingApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    // Upload files
    upload: builder.mutation({
      query: (body) => ({
        url: `upload/${crypto.randomUUID()}`,
        method: 'POST',
        body,
      }),
    }),

    // Submit refitting job
    submitRefitting: builder.mutation({
      query: (body) => ({
        url: `submitRefitting/${body.id}`,
        method: 'POST',
        body,
      }),
    }),

    // Check job status
    refittingStatus: builder.query({
      query: (jobId) => ({
        url: `refreshRefitting/${jobId}`,
      }),
    }),

    // Download results
    downloadResults: builder.query({
      query: ({ jobId, filename }) => ({
        url: `refitting/download/${jobId}/${filename}`,
      }),
    }),

    // Direct file access
    status: builder.query({
      query: (id) => ({
        url: `data/output/${id}/status.json`,
      }),
    }),

    // Multi-job status refresh
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
            const { status, params } = job;
            console.log("params --- ", params);
            return {
              jobName: params.jobName || 'Refitting Job',
              ...status,
            };
          });
      },
    }),
  }),
});

export const {
  useUploadMutation,
  useSubmitRefittingMutation,
  useRefittingStatusQuery,
  useDownloadResultsQuery,
  useStatusQuery,
  useMultiJobStatusQuery,
} = inputFormApiSlice;
