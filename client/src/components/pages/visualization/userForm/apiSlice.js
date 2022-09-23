import { visualizationApiSlice } from '../../../../services/store/rootApi';

export const userFormApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    visualizationUserUpload: builder.mutation({
      query: (formData) => ({
        url: 'upload',
        method: 'POST',
        body: formData,
      }),
    }),
    submitQueue: builder.mutation({
      query: (data) => ({
        url: 'queue',
        method: 'POST',
        body: data,
      }),
    }),
    profilerExtraction: builder.mutation({
      query: (data) => ({
        url: 'profilerExtraction',
        method: 'POST',
        body: data,
      }),
      transformResponse: (data) => {
        const matrixList = data.matrixList.map(
          ({ Profile_Type, Matrix_Size, ...e }) => ({
            ...e,
            profile: Profile_Type,
            matrix: Matrix_Size,
          })
        );
        return { ...data, matrixList };
      },
    }),
    userMatrix: builder.mutation({
      query: (params) => ({
        url: 'mutational_spectrum_options',
        params,
      }),
    }),
  }),
});

export const {
  useVisualizationUserUploadMutation,
  useSubmitQueueMutation,
  useProfilerExtractionMutation,
  useUserMatrixMutation,
} = userFormApiSlice;
