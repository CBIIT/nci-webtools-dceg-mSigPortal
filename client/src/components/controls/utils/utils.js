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
  '#BA6338',
  '#3C5488',
  '#7A65A5',
  '#F39B7F',
  '#80796B',
  '#00A1D5',
  '#6A6599',
  '#DC0000',
  '#8491B4',
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
