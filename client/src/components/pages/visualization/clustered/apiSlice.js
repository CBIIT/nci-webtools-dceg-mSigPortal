import { visualizationApiSlice } from '../../../../services/store/rootApi';
import Rainfall from '../../../controls/plotly/clusteredIdentification/rainfall';

export const clusteredApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    clustered: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        try {
          const { data } = await fetchWithBQ(
            '/cluster?' + new URLSearchParams(_arg)
          );
          const genome = data[0].genome;
          const { data: refgenome } = await fetchWithBQ(
            '/refgenome?' + new URLSearchParams({ genome })
          );
          const genomeInfo = refgenome.reduce((obj, e) => {
            const { chr, ...rest } = e;
            obj[chr] = rest;
            return obj;
          }, {});

          const plotAndTable = Rainfall(data, genomeInfo);

          return {
            data: plotAndTable,
          };
        } catch (error) {
          return { error };
        }
      },
    }),
  }),
});

export const { useClusteredQuery } = clusteredApiSlice;
