import { visualizationApiSlice } from '../../../../services/store/rootApi';
import mutationalPatternBar from '../../../controls/plotly/mutationalPattern/mutationalPatternBar';
import mutationalPatternScatter from '../../../controls/plotly/mutationalPattern/mutationalPatternScatter';
import { groupBy } from 'lodash';
import { defaultMatrix2 } from '../../../../services/utils';

export const mutationalPatternApiSlice2 = visualizationApiSlice.injectEndpoints(
  {
    endpoints: (builder) => ({
      mutationalPatternScatter: builder.query({
        query: ({ proportion, pattern, ...params }) => ({
          url: 'mutational_spectrum',
          params,
        }),
        transformResponse: (data, meta, arg) => {
          return mutationalPatternScatter(data, arg);
        },
      }),
      pattern: builder.query({
        query: (params) => ({
          url: 'pattern',
          params,
        }),

        transformResponse: (data) => {
          return mutationalPatternBar(data);
        },
      }),
    }),
  }
);

export const { useMutationalPatternScatterQuery, usePatternQuery } =
  mutationalPatternApiSlice2;
