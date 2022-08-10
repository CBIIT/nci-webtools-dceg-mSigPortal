import { visualizationApiSlice } from '../../../../services/store/rootApi';
import pcBetweenSamples from '../../../controls/plotly/profileComparision/pcBetweenSamples';
import pcToReferenceSignatures from '../../../controls/plotly/profileComparision/pcToReferenceSignatures';

export const profileComparision2 = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    pcBetweenSamples: builder.query({
      query: ({ proportion, pattern, ...params }) => ({
        url: 'seqmatrix',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        return pcBetweenSamples(data, arg);
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
