import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsLandscape from '../../../controls/plotly/msLandscape/msLandscape';

export const msLandscapeApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msLandscapePlot: builder.query({
      query: ({ variableData, ...params }) => ({
        url: 'signature_landscape',
        params,
      }),
      transformResponse: (data, meta, params) => {
        if (data.output?.error || data.output?.uncaughtError) {
          throw new Error(data.output.error || data.output.uncaughtError);
        }
        const { cosineData, exposureData, dendrogram } = data.output;
        if (cosineData && exposureData) {
          return MsLandscape(
            cosineData,
            exposureData,
            params?.variableData || [],
            dendrogram
          );
        } else {
          throw new Error('No data available');
        }
      },
    }),
  }),
});

export const { useMsLandscapePlotQuery } = msLandscapeApiSlice;
