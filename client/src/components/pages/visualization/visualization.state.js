import { atom, selector } from 'recoil';
import axios from 'axios';

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
