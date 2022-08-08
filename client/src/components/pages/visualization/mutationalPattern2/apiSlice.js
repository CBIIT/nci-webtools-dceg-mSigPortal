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

          const tmpdata0 = Object.values(
            groupBy(data, (e) => `${e.study}_${e.sample}`)
          ).map((samples) => {
            return {
              study: `${samples[0].study}`,
              sample: `${samples[0].sample}`,
              mutationType: `${samples[0].mutationType}`,
              mutation: `${samples[0].mutations}`,
              profile: `${samples[0].profile}`,
              strategy: `${samples[0].strategy}`,
              total: samples.reduce((acc, e) => acc + e.mutations, 0),
            };
          });
          console.log('TMP DATA 0 ----');
          console.log(tmpdata0);
          // const tmpdata0flat = tmpdata0.flat();
          // console.log(tmpdata0flat);

          const mutationTypeFilter = data.filter(
            (e) => e.mutationType.substring(2, 5) === type
          );
          console.log(mutationTypeFilter);
          const groupByStudySampleType = groupBy(
            mutationTypeFilter,
            (e) => `${e.study}_${e.sample}_${e.mutationType}`
          );
          console.log(groupByStudySampleType);
          const tmpdata1 = Object.values(groupByStudySampleType).map(
            (samples) => {
              return {
                study: `${samples[0].study}`,
                sample: `${samples[0].sample}`,
                mutationType: `${samples[0].mutationType}`,
                mutation: `${samples[0].mutations}`,
                profile: `${samples[0].profile}`,
                strategy: `${samples[0].strategy}`,
                n0: samples.reduce((acc, e) => acc + e.mutations, 0),
              };
            }
          );
          console.log('TMP DATA 1 ----');
          console.log(tmpdata1);
          // const tmpdata1flat = tmpdata1.flat();
          // console.log(tmpdata1flat);

          const mutationTypeSubTypesFilter = data.filter((e) =>
            subtype1 === 'N'
              ? e.mutationType.substring(2, 5) === type &&
                e.mutationType.substring(6, 7) === subtype2
              : e.mutationType.substring(2, 5) === type &&
                e.mutationType.substring(0, 1) === subtype1
          );
          console.log(mutationTypeSubTypesFilter);
          const groupByStudySampleTypeFilter = groupBy(
            mutationTypeSubTypesFilter,
            (e) => `${e.study}_${e.sample}_${e.mutationType}`
          );
          console.log(groupByStudySampleTypeFilter);

          const tmpdata2 = Object.values(groupByStudySampleTypeFilter).map(
            (samples) => {
              return {
                study: `${samples[0].study}`,
                sample: `${samples[0].sample}`,
                mutationType: `${samples[0].mutationType}`,
                mutation: `${samples[0].mutations}`,
                profile: `${samples[0].profile}`,
                strategy: `${samples[0].strategy}`,
                n1: samples.reduce((acc, e) => acc + e.mutations, 0),
              };
            }
          );
          console.log('TMP DATA 2 ----');
          console.log(tmpdata2);
          // const tmpdata2flat = tmpdata2.flat();
          // console.log(tmpdata2flat);

          function merge(a, b, prop) {
            var reduced = a.filter(
              (aitem) => !b.find((bitem) => aitem[prop] === bitem[prop])
            );
            return reduced.concat(b);
          }
          const result1 = merge(tmpdata0, tmpdata1, 'sample');
          console.log(result1);
          const result2 = merge(result1, tmpdata2, 'sample');
          console.log(result2);

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
