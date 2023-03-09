import { visualizationApiSlice } from '../../../../services/store/rootApi';
import SBS6 from '../../../controls/plotly/mutationalProfiles/sbs6';
import SBS24 from '../../../controls/plotly/mutationalProfiles/sbs24';
import SBS96 from '../../../controls/plotly/mutationalProfiles/sbs96';
import SBS192 from '../../../controls/plotly/mutationalProfiles/sbs192';
import SBS288 from '../../../controls/plotly/mutationalProfiles/sbs288';
import SBS384 from '../../../controls/plotly/mutationalProfiles/sbs384';
import SBS1536 from '../../../controls/plotly/mutationalProfiles/sbs1536';
import DBS78 from '../../../controls/plotly/mutationalProfiles/dbs78';
import DBS186 from '../../../controls/plotly/mutationalProfiles/dbs186';
import ID83 from '../../../controls/plotly/mutationalProfiles/id83';
import ID28 from '../../../controls/plotly/mutationalProfiles/id28';
import ID415 from '../../../controls/plotly/mutationalProfiles/id415';

export const mutationalProfilesApiSlice = visualizationApiSlice.injectEndpoints(
  {
    endpoints: (builder) => ({
      mutationalProfiles: builder.query({
        query: (params) => ({ url: 'mutational_spectrum', params }),
        transformResponse: (data, meta, args) => {
          const { profile, matrix, sample } = args;
          const profileMatrix = profile + matrix;

          if (profileMatrix == 'SBS6') {
            return SBS6(data);
          } else if (profileMatrix == 'SBS24') {
            return SBS24(data);
          } else if (profileMatrix == 'SBS96') {
            return SBS96(data);
          } else if (profileMatrix == 'SBS192') {
            return SBS192(data);
          } else if (profileMatrix == 'SBS288') {
            return SBS288(data);
          } else if (profileMatrix == 'SBS384') {
            return SBS384(data);
          } else if (profileMatrix == 'SBS1536') {
            return SBS1536(data);
          } else if (profileMatrix == 'DBS78') {
            return DBS78(data);
          } else if (profileMatrix == 'DBS186') {
            return DBS186(data);
          } else if (profileMatrix == 'ID83') {
            return ID83(data);
          } else if (profileMatrix == 'ID28') {
            return ID28(data);
          } else if (profileMatrix == 'ID415') {
            return ID415(data);
          } else {
            throw `Unsupported profile and matrix: ${profile} - ${matrix}`;
          }
        },
      }),
    }),
  }
);

export const { useMutationalProfilesQuery } = mutationalProfilesApiSlice;
