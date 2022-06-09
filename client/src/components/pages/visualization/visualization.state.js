import { atom, selector } from 'recoil';
import axios from 'axios';

export const defaultPublicFormState = {
  study: { label: 'PCAWG', value: 'PCAWG' },
  cancer: { label: 'Lung-AdenoCA', value: 'Lung-AdenoCA' },
  strategy: { label: 'WGS', value: 'WGS' },
};

export const publicFormState = atom({
  key: 'visualization.publicForm.state',
  default: defaultPublicFormState,
});

export const publicDataSelector = selector({
  key: 'visualization.publicData',
  get: async ({ get }) => {
    const { data } = await axios.post('web/getFileS3', {
      path: 'Others/json/Visualization-Public.json',
    });

    return data;
  },
});

export const publicFormOptions = selector({
  key: 'visualization.form.publicOptions',
  get: ({ get }) => {
    const data = get(publicDataSelector);
    const { study, cancer } = get(publicFormState);

    const studyOptions = [...new Set(data.map((d) => d.Study))].map((e) => ({
      label: e,
      value: e,
    }));

    const cancerOptions = [
      ...new Set(
        data.filter((d) => d.Study == study.value).map((d) => d.Cancer_Type)
      ),
    ].map((e) => ({ label: e, value: e }));

    const strategyOptions = [
      ...new Set(
        data
          .filter(
            (d) => d.Study == study.value && d.Cancer_Type == cancer.value
          )
          .map((d) => d.Dataset)
      ),
    ].map((e) => ({ label: e, value: e }));

    return {
      studyOptions,
      cancerOptions,
      strategyOptions,
    };
  },
});

export const defaultState = {
  source: 'public',
  inputFormat: 'vcf',
  selectedGenome: 'GRCh37',
  experimentalStrategy: 'WGS',
  mutationSplit: 'False',
  collapseSample: 'False',
  mutationFilter: '',
  queueMode: false,
  email: '',
  study: '',
  studyOptions: [],
  cancerType: '',
  cancerTypeOptions: [],
  pubExperimentalStrategy: '',
  pubExperimentOptions: [],
  pDataOptions: [],
  openSidebar: true,
  storeFilename: '',
  bedFilename: '',
  submitted: false,
  exampleData: 'assets/exampleInput/demo_input_multi.vcf.gz',
  bedData: 'assets/exampleInput/demo_input_bed.bed',
  loadingPublic: false,
  loading: {
    active: false,
    content: null,
    showIndicator: false,
  },
  queueExpired: false,
  error: '',
  projectID: '',
  displayTab: 'instructions',
  downloads: {},
  svgList: {},
  matrixList: {},
  statistics: '',
  profileOptions: [],
};

export const visualizationState = atom({
  key: 'visualization.state',
  default: defaultState,
});
// export const visualizationState = selector({
//   key: 'visualization.state',
//   get: async ({ get }) => {
//     try {
//       // const { data } = await axios.post('web/visualizationWrapper', {
//       //   fn: 'getTreeLeaf',
//       // });

//       // return data.output;
//       return defaultState;
//     } catch (error) {
//       return defaultState;
//     }
//   },
// });
