import { catalogApiSlice } from '../../../../../services/store/rootApi';
import { groupBy } from 'lodash';
import SBS96 from '../../../../controls/plotly/mutationalProfiles/sbs96';
import SBS192 from '../../../../controls/plotly/mutationalProfiles/sbs192';
import SBS288 from '../../../../controls/plotly/rsProfile/sbs288';
import SBS1536 from '../../../../controls/plotly/rsProfile/sbs1536';
import DBS78 from '../../../../controls/plotly/rsProfile/dbs78';
import ID83 from '../../../../controls/plotly/rsProfile/id83';

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
        console.log(signatureName);
        console.log(data);
        if (profileMatrix === 'SBS96') {
          return SBS96(data, signatureName, 'rsProfile');
        } else if (profileMatrix === 'SBS192') {
          return SBS192(data, signatureName, 'rsProfile');
        } else if (profileMatrix === 'SBS288') {
          return SBS288(data, signatureName);
        } else if (profileMatrix === 'SBS1536') {
          return SBS1536(data, signatureName);
        } else if (profileMatrix === 'ID83') {
          return ID83(data, signatureName);
        } else if (profileMatrix === 'DBS78') {
          return DBS78(data, signatureName);
        }
      },
    }),
  }),
});

export const { useRsProfileOptionsQuery, useRsProfileDataQuery } =
  rsProfileApiSlice;
