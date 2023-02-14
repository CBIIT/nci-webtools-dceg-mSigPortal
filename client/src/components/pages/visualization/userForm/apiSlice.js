import { visualizationApiSlice } from '../../../../services/store/rootApi';

export const userFormApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    visualizationUserUpload: builder.mutation({
      query: (formData) => ({
        url: `upload/${crypto.randomUUID()}`,
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
    submitVisualization: builder.mutation({
      query: (data) => ({
        url: 'submitVisualization',
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
    exampleHeader: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        try {
          const data = await (
            await fetch(
              `https://raw.githubusercontent.com/xtmgah/mSigPortal/master/Input_Format/Header_${_arg}.txt`
            )
          ).text();
          return { data };
        } catch (error) {
          console.log(error);
          return { error };
        }
      },
    }),
  }),
});

export const {
  useVisualizationUserUploadMutation,
  useSubmitQueueMutation,
  useSubmitVisualizationMutation,
  useUserMatrixMutation,
  useExampleHeaderQuery,
} = userFormApiSlice;
