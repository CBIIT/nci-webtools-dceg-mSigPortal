import * as d3 from 'd3';
export const customStyles = {
  control: (base, state) => ({
    ...base,
    background: '#f1e4ef',
    // match with the menu
    borderRadius: state.isFocused ? '3px 3px 0 0' : 3,
    // Overwrittes the different states of border
    borderColor: state.isFocused ? '#f1e4ef' : '#8e4b86',
    // Removes weird border around container
    boxShadow: state.isFocused ? null : null,
    '&:hover': {
      // Overwrittes the different states of border
      borderColor: state.isFocused ? '#8e4b86' : '#f1e4ef',
    },
  }),
  menu: (base) => ({
    ...base,
    // override border radius to match the box
    borderRadius: 0,
    // kill the gap
    marginTop: 0,
  }),
  menuList: (base) => ({
    ...base,
    // kill the white space on first and last option
    padding: 0,
  }),
};

export function isNumber(n) {
  return !isNaN(parseFloat(n)) && !isNaN(n - 0);
}

export const readFile = (file) => {
  let fileReader = new FileReader();
  return new Promise((resolve, reject) => {
    fileReader.onload = () => {
      resolve(fileReader.result);
    };
    fileReader.onerror = reject;
    fileReader.readAsText(file);
  });
};

export const parseMatrix = (text) => {
  let arr = text.split('\n');

  let result = [];
  for (var i = 1; i < arr.length - 1; i++) {
    let data = arr[i].split(/\t|\s+/);

    let dataObject;
    if (data.length === 3) {
      if (!isNumber(data[0].value1) && isNumber(data[1].value1)) {
        dataObject = {
          sample: data[0],
          value1: data[1],
          value2: data[2],
        };
      } else if (isNumber(data[0].value1) && !isNumber(data[1].value1)) {
        dataObject = {
          sample: data[0],
          value1: data[2],
          value2: data[1],
        };
      } else {
        dataObject = {
          sample: data[0],
          value1: data[1],
          value2: data[2],
        };
      }
    } else {
      dataObject = {
        sample: data[0],
        value1: data[1],
      };
    }
    result.push(dataObject);
  }
  return result;
};

export const checkKey = (obj, keyName) => {
  if (Object.keys(obj).indexOf(keyName) !== -1) {
    return true;
  } else {
    return false;
  }
};

export const getRandomColor = () => {
  var letters = '0123456789ABCDEF'.split('');
  var color = '#';

  for (var i = 0; i < 6; i++) {
    // color += letters[Math.floor(Math.random() * 16)];
    color += letters[Math.floor(Math.random() * 16)];
  }

  return color;
};

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

export const colorPallet = {
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
  110: '#5E1855',
  '-others': '#cececa',
};

export const mapOrder = (order, key) => (a, b) =>
  order.indexOf(a[key]) > order.indexOf(b[key]) ? 1 : -1;
