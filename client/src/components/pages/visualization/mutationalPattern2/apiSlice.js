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
          const { pattern } = arg;
          const type =
            pattern.substring(1, 2) +
            pattern.substring(3, 4) +
            pattern.substring(5, 6);
          const subtype1 = pattern.substring(0, 1);
          const subtype2 = pattern.substring(2, 3);
          const pattern1 = type + 'context';
          const pattern2 = pattern + ' other context';

          const tmpdata0 = Object.values(
            groupBy(data, (e) => `${e.study}_${e.sample}`)
          ).map((samples) => {
            return {
              study: `${samples[0].study}`,
              sample: `${samples[0].sample}`,
              type: type,
              total: samples.reduce((acc, e) => acc + e.mutations, 0),
            };
          });

          const mutationTypeFilter = data.filter(
            (e) => e.mutationType.substring(2, 5) === type
          );
          const groupByStudySampleType = groupBy(
            mutationTypeFilter,
            (e) => `${e.study}_${e.sample}`
          );
          const tmpdata1 = Object.values(groupByStudySampleType).map(
            (samples) => {
              return {
                study: `${samples[0].study}`,
                sample: `${samples[0].sample}`,
                n0: samples.reduce((acc, e) => acc + e.mutations, 0),
              };
            }
          );

          const mutationTypeSubTypesFilter = data.filter((e) =>
            subtype1 === 'N'
              ? e.mutationType.substring(2, 5) === type &&
                e.mutationType.substring(6, 7) === subtype2
              : e.mutationType.substring(2, 5) === type &&
                e.mutationType.substring(0, 1) === subtype1
          );
          const groupByStudySampleTypeFilter = groupBy(
            mutationTypeSubTypesFilter,
            (e) => `${e.study}_${e.sample}`
          );

          const tmpdata2 = Object.values(groupByStudySampleTypeFilter).map(
            (samples) => {
              return {
                study: `${samples[0].study}`,
                sample: `${samples[0].sample}`,
                type: type,
                n1: samples.reduce((acc, e) => acc + e.mutations, 0),
              };
            }
          );

          const merge = (a1, a2) =>
            a1.map((itm) => ({
              ...a2.find(
                (item) =>
                  item.sample === itm.sample && item.study === itm.study && item
              ),
              ...itm,
            }));

          const result0 = merge(tmpdata1, tmpdata0);
          const result1 = merge(tmpdata2, result0);

          const result = result1.map((e) => {
            let n2 = e.n0 - e.n1;

            return {
              study: e.study,
              sample: e.sample,
              total: e.total,
              type: pattern,
              n0: e.n0,
              n1: e.n1 / e.total,
              n2: n2 / e.total,
            };
          });
          console.log(result);

          return mutationalPatternScatter(
            result,
            type,
            subtype1,
            subtype2,
            pattern1,
            pattern2
          );
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
