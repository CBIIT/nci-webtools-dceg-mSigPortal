import { visualizationApiSlice } from '../../../../services/store/rootApi';
import Rainfall from '../../../controls/plotly/clusteredIdentification/rainfall';

export const clusteredApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    clustered: builder.query({
      query: (params) => ({
        url: 'cluster',
        params,
      }),
      transformResponse: (data) => {
        const rainfallPlot = Rainfall(
          data.filter((e) => e.sample == 'SC420396')
        );
        return { rainfall: rainfallPlot, table: data };
      },
    }),
  }),
});

export const { useClusteredQuery } = clusteredApiSlice;
