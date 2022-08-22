import { visualizationApiSlice } from '../../../../services/store/rootApi';
import pcBetweenSamples_SBS from '../../../controls/plotly/profileComparision/pcBetweenSamples_SBS';
import pcBetweenSamples_DBS from '../../../controls/plotly/profileComparision/pcBetweenSamples_DBS';
import pcBetweenSamples_ID from '../../../controls/plotly/profileComparision/pcBetweenSamples_ID';
import { groupBy } from 'lodash';

export const profilerSummaryApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    profileComparisonWithin: builder.query({
      query: (params) => ({
        url: 'seqmatrix',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        const samples = arg.sample.split(',');
        const groupBySample = groupBy(data, 'sample');
        const sample1 = groupBySample[samples[0]].flat();
        const sample2 = groupBySample[samples[1]].flat();
        if (arg.profile === 'SBS') {
          return pcBetweenSamples_SBS(samples, sample1, sample2, 'samples');
        } else if (arg.profile === 'DBS') {
          return pcBetweenSamples_DBS(samples, sample1, sample2, 'samples');
        } else {
          return pcBetweenSamples_ID(samples, sample1, sample2, 'samples');
        }
      },
    }),
    profileComparisonReference1: builder.query({
      query: (params) => ({
        url: 'seqmatrix',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        console.log(data);
        console.log(arg);
        return { ...data, arg };
      },
    }),
    profileComparisonReference2: builder.query({
      query: (params) => ({
        url: 'signature',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        console.log(data);
        console.log(arg);
        return { ...data, arg };
      },
    }),
    profileComparisonPublic: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    pcSignatureSets: builder.query({
      query: (params) => ({ url: 'signature', params }),
      transformResponse: (data) =>
        [...new Set(data.map((e) => e.signatureSetName))]
          .sort((a, b) =>
            a.localeCompare(b, undefined, { sensitivity: 'base' })
          )
          .map((e) => ({ label: e, value: e })),
    }),
    pcSignatureNames: builder.query({
      query: (params) => ({ url: 'signature', params }),
      transformResponse: (data) =>
        [...new Set(data.map((e) => e.signatureName))].sort((a, b) =>
          a.localeCompare(b, undefined, { sensitivity: 'base' })
        ),
    }),
  }),
});

export const {
  useProfileComparisonWithinQuery,
  useProfileComparisonReference1Query,
  useProfileComparisonReference2Query,
  useProfileComparisonPublicQuery,
  usePcSignatureSetsQuery,
  usePcSignatureNamesQuery,
} = profilerSummaryApiSlice;
