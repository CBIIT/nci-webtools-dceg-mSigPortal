import { catalogApiSlice } from '../../../../../services/store/rootApi';
import sbs96 from '../../../../controls/plotly/profileComparison/sbs96';
import sbs192 from '../../../../controls/plotly/profileComparison/sbs192';
import dbs78 from '../../../../controls/plotly/profileComparison/dbs78';
import id83 from '../../../../controls/plotly/profileComparison/id83';

export const rsComparisonApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    rsComparison: builder.query({
      query: (params) => ({
        url: 'mutational_signature',
        params,
      }),
      transformResponse: (data, meta, args) => {
        const signatureSetNames = args.signatureSetName.split(';');
        const signatureNames = args.signatureName.split(';');
        const data1 = data.filter(
          (e) =>
            e.signatureSetName == signatureSetNames[0] &&
            e.signatureName == signatureNames[0]
        );
        const data2 = data.filter(
          (e) =>
            e.signatureSetName == signatureSetNames[1] &&
            e.signatureName == signatureNames[1]
        );

        console.log('--RS Comparison:');
        console.log(data1);
        console.log(data2);
        if (args.profile == 'SBS' && args.matrix == '96')
          return sbs96(data1, data2);
        else if (args.profile == 'SBS' && args.matrix == '192')
          return sbs192(data1, data2);
        else if (args.profile == 'DBS') return dbs78(data1, data2);
        else if (args.profile == 'ID') return id83(data1, data2);
        else throw Error(`Profile ${args.profile} is not supported`);
      },
    }),
  }),
});

export const { useRsComparisonQuery } = rsComparisonApiSlice;
