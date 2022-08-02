import { visualizationApiSlice } from '../../../../services/store/rootApi';
import mutationalPattern from '../../../controls/plotly/mutationalPattern/mutationalPattern';
import { groupBy } from 'lodash';

export const mutationalPatternApiSlice2 = visualizationApiSlice.injectEndpoints(
  {
    endpoints: (builder) => ({
      mutationalPattern2: builder.query({
        query: ({ proportion, pattern, ...params }) => ({
          url: 'seqmatrix',
          params,
        }),
        transformResponse: (data, meta, arg) => {
          console.log(data);
          const { proportion, pattern } = arg;
          console.log(proportion);
          console.log(pattern);
          const groupByProfile = groupBy(data, 'profile');

          return mutationalPattern(groupByProfile);
        },
      }),
    }),
  }
);

export const { useMutationalPattern2Query } = mutationalPatternApiSlice2;
