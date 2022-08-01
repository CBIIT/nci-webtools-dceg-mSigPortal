import { visualizationApiSlice } from '../../../../services/store/rootApi';
import mutationalPattern from '../../../controls/plotly/mutationalPattern/mutationalPattern';
import { groupBy } from 'lodash';

export const mutationalPatternApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    mutationalPattern: builder.query({
      query: (params) => ({
        url: 'seqmatrix',
        body: params,
      }),
      transformResponse: (data) => {
        console.log(data);
        const groupByProfile = groupBy(data, 'profile');

        return mutationalPattern(groupByProfile);
      },
    }),
  }),
});

export const { useMutationalPatternQuery } = mutationalPatternApiSlice;
