import { visualizationApiSlice } from '../../../../services/store/rootApi';
import mutationalPattern from '../../../controls/plotly/mutationalPattern/mutationalPattern';
import { groupBy } from 'lodash';

export const mutationalPatternApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    mutationalPattern: builder.query({
      query: ({ proportion, pattern, ...params }) => ({
        url: 'seqmatrix',
        body: params,
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
});

export const { useMutationalPatternQuery } = mutationalPatternApiSlice;
