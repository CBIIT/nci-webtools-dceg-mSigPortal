import { apiSlice } from '../../../../services/apiSlice';

export const vissualizationApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    userFormUpload: builder.mutation({
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

    submitWebUser: builder.mutation({
      query: (data) => ({
        url: 'profilerExtraction',
        method: 'POST',
        body: data,
      }),
      transformResponse: (data) => {
        const svgList = data.svgList.map(
          ({ Sample_Name, Profile_Type, Matrix_Size, ...e }) => ({
            ...e,
            sample: Sample_Name,
            profileType: Profile_Type,
            matrixSize: Matrix_Size,
          })
        );
        const matrixList = data.matrixList.map(
          ({ Profile_Type, Matrix_Size, ...e }) => ({
            ...e,
            profileType: Profile_Type,
            matrixSize: Matrix_Size,
          })
        );
        return { ...data, svgList, matrixList };
      },
    }),
  }),
});

export const {
  useUserFormUploadMutation,
  useSubmitQueueMutation,
  useSubmitWebUserMutation,
} = vissualizationApiSlice;
