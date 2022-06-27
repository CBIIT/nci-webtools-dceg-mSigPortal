import { apiSlice } from '../../../../services/apiSlice';
import { groupBy } from 'lodash';

export const vissualizationApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    submitWebPublic: builder.mutation({
      query: (data) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: { fn: 'getPublicData', args: data },
      }),
      transformResponse: (data) => {
        const svgList = data.svgList.map(({ Sample, ...e }) => ({
          ...e,
          // sample: Sample,
          profileType: e.Profile.match(/[a-z]+/gi)[0],
          matrixSize: e.Profile.match(/\d+/gi)[0],
        }));

        return { ...data, svgList };
      },
    }),

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
  }),
});

export const { useSubmitWebPublicMutation, useGetPublicDataOptionsQuery } =
  vissualizationApiSlice;
