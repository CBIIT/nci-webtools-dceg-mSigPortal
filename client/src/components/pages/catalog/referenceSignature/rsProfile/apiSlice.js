import { catalogApiSlice } from '../../../../../services/store/rootApi';
import { groupBy } from 'lodash';
import SBS96 from '../../../../controls/plotly/mutationalProfiles/sbs96';
import SBS192 from '../../../../controls/plotly/mutationalProfiles/sbs192';

export const rsProfileApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    rsProfileOptions: builder.query({
      query: (params) => ({
        url: 'mutational_signature_options',
        params,
      }),
    }),
    rsProfileData: builder.query({
      query: (params) => ({
        url: 'mutational_signature',
        params,
      }),
      transformResponse: (data, meta, args) => {
        console.log(args);
        const { profile, matrix, signatureName } = args;
        const profileMatrix = profile + matrix;
        console.log(profileMatrix);
        console.log(data);
        if (profileMatrix === 'SBS96') {
          return SBS96(data, signatureName, 'rsProfile');
        } else if (profileMatrix === 'SBS192') {
          return SBS192(data, signatureName, 'rsProfile');
        }
        return { data };
      },
    }),
  }),
});

export const { useRsProfileOptionsQuery, useRsProfileDataQuery } =
  rsProfileApiSlice;
