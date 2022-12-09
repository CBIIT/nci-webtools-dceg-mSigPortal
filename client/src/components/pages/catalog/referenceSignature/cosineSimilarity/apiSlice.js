import { catalogApiSlice } from '../../../../../services/store/rootApi';
import CosineSimilarity from '../../../../controls/plotly/cosineSimilarity/cosineSimilarity';

export const cosineSimilarityApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    cosineSimilarity: builder.query({
      query: (params) => ({
        url: 'signature_cosine_similarity',
        params,
      }),
      transformResponse: (data) => {
        const { output } = data;
        return { ...CosineSimilarity(output), original: output.original };
      },
    }),
  }),
});

export const { useCosineSimilarityQuery } = cosineSimilarityApiSlice;
