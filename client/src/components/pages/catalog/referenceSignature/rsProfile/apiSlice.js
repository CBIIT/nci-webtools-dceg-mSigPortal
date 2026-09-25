import { catalogApiSlice } from '@/services/store/rootApi';
import SBS96 from '@/components/controls/plotly/rsProfile/sbs96';
import SBS192 from '@/components/controls/plotly/rsProfile/sbs192';
import SBS288 from '@/components/controls/plotly/rsProfile/sbs288';
import SBS1536 from '@/components/controls/plotly/rsProfile/sbs1536';
import DBS78 from '@/components/controls/plotly/rsProfile/dbs78';
import ID83 from '@/components/controls/plotly/rsProfile/id83';
import RS32 from '@/components/controls/plotly/rsProfile/rs32';
import CN48 from '@/components/controls/plotly/rsProfile/cn48';

export const rsProfileApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    rsProfileOptions: builder.query({
      query: (params) => ({
        url: 'mutational_signature_options',
        params,
      }),
    }),

    rsProfilePlot: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        const params = _arg.params;
        try {
          const res = await Promise.all(
            params.map((e) =>
              fetchWithBQ('/mutational_signature?' + new URLSearchParams(e))
            )
          );

          const plotData = res.map((e, i) => {
            if (params[i].profile + params[i].matrix === 'SBS96') {
              return SBS96(res[i].data, params[i].signatureName, 'rsProfile');
            } else if (params[i].profile + params[i].matrix === 'SBS192') {
              return SBS192(res[i].data, params[i].signatureName, 'rsProfile');
            } else if (params[i].profile + params[i].matrix === 'SBS288') {
              return SBS288(res[i].data, params[i].signatureName);
            } else if (params[i].profile + params[i].matrix === 'SBS1536') {
              return SBS1536(res[i].data, params[i].signatureName);
            } else if (params[i].profile + params[i].matrix === 'ID83') {
              return ID83(res[i].data, params[i].signatureName);
            } else if (params[i].profile + params[i].matrix === 'DBS78') {
              return DBS78(res[i].data, params[i].signatureName);
            } else if (params[i].profile + params[i].matrix === 'RS32') {
              return RS32(res[i].data, params[i].signatureName);
            } else if (params[i].profile + params[i].matrix === 'CN48') {
              return CN48(res[i].data, params[i].signatureName);
            }
          });

          if (plotData.length) {
            return { data: plotData };
          } else {
            return { data: [] };
          }
        } catch (error) {
          return { error };
        }
      },
    }),
  }),
});

export const { useRsProfileOptionsQuery, useRsProfilePlotQuery } =
  rsProfileApiSlice;
