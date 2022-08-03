import { visualizationApiSlice } from '../../../../services/store/rootApi';
import mutationalPatternBar from '../../../controls/plotly/mutationalPattern/mutationalPatternBar';
import mutationalPatternScatter from '../../../controls/plotly/mutationalPattern/mutationalPatternScatter';
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
          const { proportion, pattern, profile, matrix } = arg;
          const regexMap = {
            SBS24: /^.{2}(.*)$/,
            SBS96: /\[(.*)\]/,
            DBS78: /^(.{2})/,
            ID83: /^(.{7})/,
          };
          const profileMatrix = profile + matrix;
          const groupByMutationType = groupBy(data, 'mutationType');
          console.log(groupByMutationType);
          const groupByMutation = data.reduce((acc, e, i) => {
            const mutationRegex = regexMap[profileMatrix];
            const mutation = e.mutationType.match(mutationRegex)[1];

            acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
            return acc;
          }, {});

          const transform = Object.entries(groupByMutation).map(
            ([mutation, data]) => ({
              mutation,
              data,
            })
          );

          console.log(transform);

          return mutationalPatternScatter(transform, proportion, pattern);
        },
      }),
    }),
  }
);

export const { useMutationalPattern2Query } = mutationalPatternApiSlice2;
