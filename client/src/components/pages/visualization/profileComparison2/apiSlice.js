import { visualizationApiSlice } from '../../../../services/store/rootApi';
import pcBetweenSamples from '../../../controls/plotly/profileComparision/pcBetweenSamples';
import pcToReferenceSignatures from '../../../controls/plotly/profileComparision/pcToReferenceSignatures';

export const profileComparision2 = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    pcBetweenSamples: builder.query({
      query: ({ params }) => ({
        url: 'seqmatrix',
        params,
      }),

      transformResponse: (data, meta, args) => {
        const { profile, sample1, sample2 } = args;
        console.log(data);
        console.log(profile);
        console.log(sample1);
        console.log(sample2);
        return pcBetweenSamples(data, profile, sample1, sample2);
      },
    }),
    pcToReferenceSignatures: builder.query({
      query: (params) => ({
        url: 'seqmatrix',
        params,
      }),

      transformResponse: (data, meta, arg) => {
        return pcToReferenceSignatures(data, arg);
      },
    }),
  }),
});

export const { usePcBetweenSamplesQuery, usePcToReferenceSignaturesQuery } =
  profileComparision2;
