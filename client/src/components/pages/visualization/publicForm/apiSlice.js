import { apiSlice } from '../../../../services/apiSlice';
import { groupBy } from 'lodash';

export const publicFormApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    getPublicDataOptions: builder.query({
      query: (_) => ({
        url: 'getFileS3',
        method: 'POST',
        body: { path: 'Others/json/Visualization-Public.json' },
      }),
      transformResponse: (data) => {
        const groupByStudy = groupBy(data, 'Study');
        const groupByCancer = Object.fromEntries(
          Object.entries(groupByStudy).map(([study, data]) => [
            study,
            Object.entries(groupBy(data, 'Cancer_Type')).reduce(
              (a, [cancer, e]) => ((a[cancer] = { Datset: e[0].Dataset }), a),
              {}
            ),
          ])
        );

        return groupByCancer;
      },
    }),
    publicSamples: builder.mutation({
      query: (params) => ({
        url: 'visualizationData',
        params,
      }),
    }),
  }),
});

export const { useGetPublicDataOptionsQuery, usePublicSamplesMutation } =
  publicFormApiSlice;
