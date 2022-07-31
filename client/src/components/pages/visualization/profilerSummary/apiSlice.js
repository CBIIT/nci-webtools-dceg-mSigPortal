import { visualizationApiSlice } from '../../../../services/store/rootApi';
import profilerSummary from '../../../controls/plotly/profilerSummary/profilerSummary';
import { groupBy } from 'lodash';

export const profilerSummaryApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    profilerSummary: builder.query({
      query: (params) => ({
        url: 'seqmatrixSummary',
        params,
      }),
      transformResponse: (data) => {
        const groupByProfile = groupBy(data, 'profile');
        const meanMutationsPerProfile = Object.values(groupByProfile)
          .map((samples) => {
            return {
              name: `${samples[0].profile}: ${samples[0].matrix}`,
              samples: samples.sort(
                (a, b) => a.logTotalMutations - b.logTotalMutations
              ),
              mean:
                samples.reduce(
                  (acc, e) => acc + parseFloat(e.meanTotalMutations),
                  0
                ) / samples.length,
            };
          })
          .sort((a, b) => b.mean - a.mean);

        return profilerSummary(meanMutationsPerProfile);
      },
    }),
  }),
});

export const { useProfilerSummaryQuery } = profilerSummaryApiSlice;
