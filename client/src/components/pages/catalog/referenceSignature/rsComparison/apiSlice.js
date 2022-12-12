import { catalogApiSlice } from '../../../../../services/store/rootApi';
import profileComparision_SBS from '../../../../controls/plotly/profileComparision/pcBetweenSamples_SBS';

export const rsComparisonApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    rsComparison: builder.query({
      query: (params) => ({
        url: 'mutational_signature',
        params,
      }),
      transformResponse: (data) => {
        // return profileComparision_SBS(data);
        return data;
      },
    }),
  }),
});

export const { useRsComparisonQuery } = rsComparisonApiSlice;
