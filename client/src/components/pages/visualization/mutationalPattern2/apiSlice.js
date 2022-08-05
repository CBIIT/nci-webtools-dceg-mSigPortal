import { visualizationApiSlice } from '../../../../services/store/rootApi';
import mutationalPatternBar from '../../../controls/plotly/mutationalPattern/mutationalPatternBar';
import mutationalPatternScatter from '../../../controls/plotly/mutationalPattern/mutationalPatternScatter';
import { groupBy } from 'lodash';

export const mutationalPatternApiSlice2 = visualizationApiSlice.injectEndpoints(
  {
    endpoints: (builder) => ({
      mutationalPatternScatter: builder.query({
        query: ({ proportion, pattern, ...params }) => ({
          url: 'seqmatrix',
          params,
        }),
        transformResponse: (data, meta, arg) => {
          console.log(data);
          const { proportion, pattern, profile, matrix } = arg;
          const type =
            pattern.substring(1, 2) +
            pattern.substring(3, 4) +
            pattern.substring(5, 6);
          const subtype1 = pattern.substring(0, 1);
          const subtype2 = pattern.substring(2, 3);

          // const regexMap = {
          //   SBS24: /^.{2}(.*)$/,
          //   SBS96: /\[(.*)\]/,
          //   DBS78: /^(.{2})/,
          //   ID83: /^(.{7})/,
          // };
          // const profileMatrix = profile + matrix;
          // const groupByMutationType = groupBy(data, 'mutationType');
          // console.log(groupByMutationType);
          // const groupBySample = groupBy(data, 'sample');
          // console.log(groupBySample);
          // const groupByMutation = data.reduce((acc, e, i) => {
          //   const mutationRegex = regexMap[profileMatrix];
          //   const mutation = e.mutationType.match(mutationRegex)[1];
          //   acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
          //   return acc;
          // }, {});

          // const transform = Object.entries(groupByMutation).map(
          //   ([mutation, data]) => {
          //     const mutationSum = data.reduce((sum, e) => e.mutations + sum, 0);
          //     return { mutation, data, mutationSum };
          //   }
          // );

          // console.log(transform);
          // const filter = transform.filter((e) => {
          //   return e.mutation === type;
          // });

          const groupByStudySample = groupBy(
            data,
            (e) => `${e.study}_${e.sample}`
          );
          console.log(groupByStudySample);
          const total = Object.values(groupByStudySample).map((samples) => {
            return {
              name: `${samples[0].study}_${samples[0].sample}`,
              total: samples.reduce((acc, e) => acc + e.mutations, 0),
            };
          });
          console.log(total);
          const mutationTypeFilter = data.filter(
            (e) => e.mutationType.substring(2, 5) === type
          );
          console.log(mutationTypeFilter);
          const groupByStudySampleType = groupBy(
            mutationTypeFilter,
            (e) => `${e.study}_${e.sample}_${e.mutationType}`
          );
          console.log(groupByStudySampleType);
          const n0 = Object.values(groupByStudySampleType).map((samples) => {
            return {
              name: `${samples[0].study}_${samples[0].sample}_${samples[0].mutationType}`,
              total: samples.reduce((acc, e) => acc + e.mutations, 0),
            };
          });
          console.log(n0);

          return mutationalPatternScatter(type, subtype1, subtype2);
        },
      }),
      pattern: builder.query({
        query: ({ proportion, ...params }) => ({
          url: 'pattern',
          params,
        }),
        transformResponse: (data, meta, arg) => {
          console.log(data);
          const { proportion } = arg;
          const transform = data.filter((e) => e.n1 >= proportion);
          console.log(transform);
          const groupByPattern = groupBy(transform, 'pattern');
          console.log(groupByPattern);
          const transform2 = Object.entries(groupByPattern).map(
            ([pattern, data]) => {
              return { pattern, data };
            }
          );
          return mutationalPatternBar(transform2);
        },
      }),
    }),
  }
);

export const { useMutationalPatternScatterQuery, usePatternQuery } =
  mutationalPatternApiSlice2;
