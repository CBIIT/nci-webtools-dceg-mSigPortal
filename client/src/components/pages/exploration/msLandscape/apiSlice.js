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
        const { cosineData, exposureData, dendrogram } = data.output;
        if (cosineData && exposureData) {
          return MsLandscape(
            cosineData,
            exposureData,
            params?.variableData || [],
            dendrogram
          );
        } else {
          return false;
        }
      },
    }),
  }),
});

export const { useMsLandscapePlotQuery } = msLandscapeApiSlice;
