export const colorPallet0 = [
  '#00A1D5',
  '#F39B7F',
  '#BA6338',
  '#8491B4',
  '#3C5488',
  '#7A65A5',
  '#80796B',
  '#6A6599',
  '#DC0000',
  '#837B8D',
  '#E7C76F',
  '#DF8F44',
  '#CDDEB7',
  '#E64B35',
  '#91D1C2',
  '#924822',
  '#B09C85',
  '#F0E685',
  '#4DBBD5',
  '#AE1F63',
  '#7E6148',
  '#612A79',
  '#5050FF',
  '#D595A7',
  '#466983',
  '#3B1B53',
  '#749B58',
  '#802268',
  '#5DB1DD',
  '#E4AF69',
  '#79AF97',
  '#00A087',
  '#374E55',
  '#B24745',
  '#C75127',
  '#CE3D32',
  '#6BD76B',
  '#D58F5C',
];

export const colorPallet1 = [
  '#1F77B4',
  '#FF7F0F',
  '#2DA02C',
  '#D62728',
  '#9467BD',
  '#8C564B',
];

const originalColorPallet = {
  1: '#4a9855',
  2: '#e2a8ab',
  3: '#40004b',
  4: '#5aa1ca',
  5: '#305d39',
  6: '#785940',
  '7a': '#6e70b7',
  '7b': '#ff7f00',
  '7c': '#fec44f',
  '7d': '#846a2a',
  8: '#cab2d6',
  9: '#f4a582',
  '10a': '#8dd3c7',
  '10b': '#5e4fa2',
  '10c': '#761429',
  11: '#9e0142',
  12: '#ffed6f',
  13: '#e41a1c',
  14: '#ffffbf',
  15: '#4d4d4d',
  16: '#513276',
  17: '#ef4c7d',
  '17a': '#df4c7d',
  '17b': '#08519c',
  18: '#b3de69',
  19: '#dfc27d',
  20: '#b2182b',
  21: '#9ecae1',
  22: '#01665e',
  23: '#d53e4f',
  24: '#1c9099',
  25: '#35978f',
  26: '#ec7014',
  27: '#f46d43',
  28: '#de77ae',
  29: '#fdae61',
  30: '#d9d9d9',
  31: '#f781bf',
  32: '#dd1c77',
  33: '#b25d7e',
  34: '#fee08b',
  35: '#fc8d59',
  36: 'yellow',
  37: '#e6f598',
  38: '#abdda4',
  39: '#636363',
  40: '#b15928',
  41: '#fccde5',
  42: '#ae017e',
  43: '#66c2a5',
  44: '#8c6bb1',
  45: '#3288bd',
  46: '#e6f598',
  47: '#bababa',
  48: '#5e4fa2',
  49: '#40004b',
  50: '#762a83',
  51: '#9970ab',
  52: '#c2a5cf',
  53: '#e7d4e8',
  54: '#fcc5c0',
  55: '#d9f0d3',
  56: '#8c510a',
  57: '#a6dba0',
  58: '#5aae61',
  59: '#1b7837',
  60: '#00441b',
  84: '#063C3C',
  85: '#AA9139',
  88: '#BB9139',
  92: '#0E1844',
  96: '#1F77B5',
  '96A': '#2F77B5',
  '96B': '#3F77B5',
  '96C': '#4F77B5',
  '96D': '#1FBBB5',
  110: '#5E1855',
  '-others': '#cececa',
};

// The new object
const additionalColorPallet = {};

// Generate random colors for each entry not in originalColorPallet
for (let i = 1; i <= 110; i++) {
  const baseKey = i.toString();

  // Add the base entry if it doesn't already exist
  if (!originalColorPallet[baseKey]) {
    additionalColorPallet[baseKey] = `#${Math.floor(
      Math.random() * 16777215
    ).toString(16)}`;
  }

  // Add the four sub-entries if they don't already exist
  for (let j = 97; j <= 100; j++) {
    const subKey = `${baseKey}${String.fromCharCode(j)}`;
    if (!originalColorPallet[subKey]) {
      additionalColorPallet[subKey] = `#${Math.floor(
        Math.random() * 16777215
      ).toString(16)}`;
    }

    const subKeyUpper = `${baseKey}${String.fromCharCode(j).toUpperCase()}`;
    if (!originalColorPallet[subKeyUpper]) {
      additionalColorPallet[subKeyUpper] = `#${Math.floor(
        Math.random() * 16777215
      ).toString(16)}`;
    }
  }
}

// Add the '-others' entry if it doesn't already exist
if (!originalColorPallet['-others']) {
  additionalColorPallet['-others'] = '#cececa';
}

// Merge the two objects using Object.assign()
export const colorPallet = Object.assign(
  {},
  originalColorPallet,
  additionalColorPallet
);

export const rs32Color = {
  'clustered_del_>10Mb': 'deeppink',
  'non-clustered_del_>10Mb': 'deeppink',
  'clustered_del_1Mb-10Mb': 'hotpink',
  'non-clustered_del_1Mb-10Mb': 'hotpink',
  'clustered_del_10-100Kb': 'lightpink',
  'non-clustered_del_10-100Kb': 'lightpink',
  'clustered_del_100Kb-1Mb': 'palevioletred',
  'non-clustered_del_100Kb-1Mb': 'palevioletred',
  'clustered_del_1-10Kb': 'lavenderblush',
  'non-clustered_del_1-10Kb': 'lavenderblush',
  'clustered_tds_>10Mb': 'saddlebrown',
  'non-clustered_tds_>10Mb': 'saddlebrown',
  'clustered_tds_1Mb-10Mb': 'sienna',
  'non-clustered_tds_1Mb-10Mb': 'sienna',
  'clustered_tds_10-100Kb': 'sandybrown',
  'non-clustered_tds_10-100Kb': 'sandybrown',
  'clustered_tds_100Kb-1Mb': 'peru',
  'non-clustered_tds_100Kb-1Mb': 'peru',
  'clustered_tds_1-10Kb': 'linen',
  'non-clustered_tds_1-10Kb': 'linen',
  'clustered_inv_>10Mb': 'rebeccapurple',
  'non-clustered_inv_>10Mb': 'rebeccapurple',
  'clustered_inv_1Mb-10Mb': 'blueviolet',
  'non-clustered_inv_1Mb-10Mb': 'blueviolet',
  'clustered_inv_10-100Kb': 'plum',
  'non-clustered_inv_10-100Kb': 'plum',
  'clustered_inv_100Kb-1Mb': 'mediumorchid',
  'non-clustered_inv_100Kb-1Mb': 'mediumorchid',
  'clustered_inv_1-10Kb': 'thistle',
  'non-clustered_inv_1-10Kb': 'thistle',
  clustered_trans: 'gray',
  'non-clustered_trans': 'gray',
  del: '#800001',
  tds: '#FF8C00',
  inv: '#6A5ACD',
  tra: '#696969',
};

export const sbsColor = {
  'C>A': '#03BCEE',
  'C>G': 'black',
  'C>T': '#E32926',
  'T>A': '#CAC9C9',
  'T>C': '#A1CE63',
  'T>G': '#EBC6C4',
};

export const id83Color = {
  '1:Del:C': { shape: '#FBBD6F', text: 'black' },
  '1:Del:T': { shape: '#FE8002', text: 'white' },
  '1:Ins:C': { shape: '#AEDD8A', text: 'black' },
  '1:Ins:T': { shape: '#35A12E', text: 'white' },
  '2:Del:R': { shape: '#FCC9B4', text: 'black' },
  '3:Del:R': { shape: '#FB8969', text: 'black' },
  '4:Del:R': { shape: '#F04432', text: 'black' },
  '5:Del:R': { shape: '#BB1A1A', text: 'white' },
  '2:Ins:R': { shape: '#CFDFF0', text: 'black' },
  '3:Ins:R': { shape: '#93C3DE', text: 'black' },
  '4:Ins:R': { shape: '#4B97C7', text: 'black' },
  '5:Ins:R': { shape: '#1863AA', text: 'white' },
  '2:Del:M': { shape: '#E1E1EE', text: 'blacl' },
  '3:Del:M': { shape: '#B5B5D6', text: 'black' },
  '4:Del:M': { shape: '#8482BC', text: 'black' },
  '5:Del:M': { shape: '#62409A', text: 'white' },
};

export const dbs78Color = {
  AC: '#09BCED',
  AT: '#0266CA',
  CC: '#9FCE62',
  CG: '#006501',
  CT: '#FF9898',
  GC: '#E22925',
  TA: '#FEB065',
  TC: '#FD8000',
  TG: '#CB98FD',
  TT: '#4C0299',
};
