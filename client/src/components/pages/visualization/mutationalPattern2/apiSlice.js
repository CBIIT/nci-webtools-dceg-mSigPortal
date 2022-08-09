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
          const pattern1 = type + ' context';
          const pattern2 = pattern + ' other context';
          console.log(type);

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
          console.log('TMPDATA 0');
          console.log(tmpdata0);
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
          console.log('TMPDATA 1');
          console.log(tmpdata1);

          function iupac(base) {
            let result = [];
            if (base === 'T' || base === 'U') {
              result = ['T'];
            } else if (base === 'M') {
              result = ['A', 'C'];
            } else if (base === 'R') {
              result = ['A', 'G'];
            } else if (base === 'W') {
              result = ['A', 'T'];
            } else if (base === 'S') {
              result = ['C', 'G'];
            } else if (base === 'Y') {
              result = ['C', 'T'];
            } else if (base === 'K') {
              result = ['G', 'T'];
            } else if (base === 'V') {
              result = ['A', 'C', 'G'];
            } else if (base === 'H') {
              result = ['A', 'C', 'T'];
            } else if (base === 'D') {
              result = ['A', 'G', 'T'];
            } else if (base === 'B') {
              result = ['C', 'G', 'T'];
            } else if (base === 'N') {
              result = ['G', 'A', 'T', 'C'];
            } else {
              result = [base];
            }
            return result;
          }
          const mutationTypeSubTypesFilter = data.filter(
            (e) =>
              e.mutationType.substring(2, 5) === type &&
              e.mutationType.substring(0, 1) === iupac(subtype1) &&
              e.mutationType.substring(6, 7) === iupac(subtype2)
          );
          // const mutationTypeSubTypesFilter = data.filter((e) =>
          //   subtype1 === 'N' && subtype2 !== 'N'
          //     ? e.mutationType.substring(2, 5) === type &&
          //       e.mutationType.substring(6, 7) === subtype2
          //     : e.mutationType.substring(2, 5) === type &&
          //       e.mutationType.substring(0, 1) === subtype1
          // );
          console.log(mutationTypeSubTypesFilter);
          const groupByStudySampleTypeFilter = groupBy(
            mutationTypeSubTypesFilter,
            (e) => `${e.study}_${e.sample}`
          );

          const tmpdata2 = Object.values(groupByStudySampleTypeFilter).map(
            (samples) => {
              return {
                study: `${samples[0].study}`,
                sample: `${samples[0].sample}`,
                cancer: `${samples[0].cancer}`,
                type: type,
                n1: samples.reduce((acc, e) => acc + e.mutations, 0),
              };
            }
          );
          console.log('TMPDATA 2');
          console.log(tmpdata2);
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

          const result = result1
            .map((e) => {
              let n2 = e.n0 - e.n1;

              return {
                study: e.study,
                sample: e.sample,
                cancer: e.cancer,
                total: e.total,
                type: pattern,
                n0: e.n0,
                n1: e.n1 / e.total,
                n2: n2 / e.total,
              };
            })
            .filter((e) => e.total > 200);
          console.log(result);

          return mutationalPatternScatter(result, type, pattern1, pattern2);
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
          const transform = data.filter((e) => e.n1 > proportion);
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
