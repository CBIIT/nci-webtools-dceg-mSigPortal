import { apiSlice } from "../../../../services/apiSlice";
import TMB from "../../../controls/plotly/tmb/tmb";

export const tmbApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    tmbPlot: builder.query({
      query: (params) => ({
        url: "explorationTmbData",
        params,
      }),
      transformResponse: (data, meta, arg) => {
        return TMB(data, arg.study);
      },
    }),
  }),
});

export const { useTmbPlotQuery } = tmbApiSlice;
