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
        return MsLandscape(
          cosineData,
          exposureData,
          params?.variableData || [],
          dendrogram
        );
      },
    }),
  }),
});

export const { useMsLandscapePlotQuery } = msLandscapeApiSlice;
