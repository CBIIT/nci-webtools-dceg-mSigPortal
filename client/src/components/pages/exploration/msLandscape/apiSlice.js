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
        console.log('ms landscape data', data);
        const { cosineData, exposureData, dendrogram } = data.output;
        console.log('cosineData', cosineData);
        console.log('exposureData', exposureData);
        console.log('dendrogram', dendrogram);
        if (cosineData && exposureData) {
          return MsLandscape(
            cosineData,
            exposureData,
            params?.variableData || [],
            dendrogram
          );
        } else {
          throw new Error(data.stdout);
        }
      },
    }),
  }),
});

export const { useMsLandscapePlotQuery } = msLandscapeApiSlice;
