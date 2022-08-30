import { visualizationApiSlice } from '../../../../services/store/rootApi';
import pcBetweenSamples_SBS from '../../../controls/plotly/profileComparision/pcBetweenSamples_SBS';
import pcBetweenSamples_DBS from '../../../controls/plotly/profileComparision/pcBetweenSamples_DBS';
import pcBetweenSamples_ID from '../../../controls/plotly/profileComparision/pcBetweenSamples_ID';
import { groupBy } from 'lodash';
import { createApi, fetchBaseQuery } from '@reduxjs/toolkit/query';

export const profilerSummaryApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    profileComparisonWithin: builder.query({
      query: (params) => ({
        url: 'mutational_spectrum',
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
    // profileComparisonReference1: builder.query({
    //   query: (params) => ({
    //     url: 'mutational_spectrum',
    //     params,
    //   }),
    //   transformResponse: (data, meta, arg) => {
    //     console.log(data);
    //     console.log(arg);
    //     return { ...data, arg };
    //   },
    // }),
    // profileComparisonReference2: builder.query({
    //   query: (params) => ({
    //     url: 'mutational_signature',
    //     params,
    //   }),
    //   transformResponse: (data, meta, arg) => {
    //     console.log(data);
    //     console.log(arg);
    //     return { ...data, arg };
    //   },
    // }),

    profileComparisonReference: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        console.log(_arg);
        // get a random user
        const { data: spectrumData, error: spectrumError } = await fetchWithBQ(
          'mutational_spectrum?' + new URLSearchParams(_arg.params_spectrum)
        );
        if (spectrumError) return { error: spectrumError };

        const { data: signatureData, error: signatureError } =
          await fetchWithBQ(
            'mutational_signature?' + new URLSearchParams(_arg.params_signature)
          );
        if (signatureError) return { error: signatureError };

        console.log(spectrumData);
        console.log(signatureData);
        // console.log(_arg);

        return { spectrumData, signatureData };
      },
      transformResponse: (data, meta, arg) => {
        console.log(data);
        console.log(meta);
        console.log(arg);
        //   const samples = [
        //   _arg.params_spectrum.sample,
        //   _arg.params_signature.signatureName,
        // ];
        // // console.log(samples);
        // if (_arg.params_spectrum.profile === 'SBS') {
        //   return pcBetweenSamples_SBS(
        //     samples,
        //     spectrumData,
        //     signatureData,
        //     'reference'
        //   );
        //   // return pcBetweenSamples_SBS({
        //   //   samples: samples,
        //   //   sample1: spectrumData,
        //   //   sample2: signatureData,
        //   //   tab: 'reference',
        //   // });
        // } else if (_arg.params_spectrum.profile === 'DBS') {
        //   return pcBetweenSamples_DBS(
        //     samples,
        //     spectrumData,
        //     signatureData,
        //     'reference'
        //   );
        // } else {
        //   return pcBetweenSamples_ID(
        //     samples,
        //     spectrumData,
        //     signatureData,
        //     'reference'
        //   );
        // }}
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
      query: (params) => ({ url: 'mutational_signature', params }),
      transformResponse: (data) =>
        [...new Set(data.map((e) => e.signatureSetName))]
          .sort((a, b) =>
            a.localeCompare(b, undefined, { sensitivity: 'base' })
          )
          .map((e) => ({ label: e, value: e })),
    }),
    pcSignatureNames: builder.query({
      query: (params) => ({ url: 'mutational_signature', params }),
      transformResponse: (data, meta, arg) =>
        [...new Set(data.map((e) => e.signatureName))].sort((a, b) =>
          a.localeCompare(b, undefined, { sensitivity: 'base' })
        ),
    }),
  }),
});

export const {
  useProfileComparisonWithinQuery,

  useProfileComparisonReferenceQuery,
  useProfileComparisonPublicQuery,
  usePcSignatureSetsQuery,
  usePcSignatureNamesQuery,
} = profilerSummaryApiSlice;
