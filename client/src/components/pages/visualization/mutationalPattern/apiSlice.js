import { visualizationApiSlice } from '../../../../services/store/rootApi';
import mutationalPatternBar from '../../../controls/plotly/mutationalPattern/mutationalPatternBar';
import mutationalPatternScatter from '../../../controls/plotly/mutationalPattern/mutationalPatternScatter';

export const mutationalPatternApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    mpeaScatter: builder.query({
      query: ({ proportion, pattern, ...params }) => ({
        url: 'mutational_spectrum',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        return mutationalPatternScatter(data, arg);
      },
    }),
    mpeaBar: builder.query({
      query: (params) => ({
        url: 'pattern',
        params,
      }),

      transformResponse: (data) => {
        if (data.length) return mutationalPatternBar(data);
        else return {};
      },
    }),
  }),
});

export const { useMpeaScatterQuery, useMpeaBarQuery } =
  mutationalPatternApiSlice;
