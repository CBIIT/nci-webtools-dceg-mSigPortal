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

export const mapOrder = (order, key) => (a, b) =>
  order.indexOf(a[key]) > order.indexOf(b[key]) ? 1 : -1;

export function groupByCustom(list, keyGetter) {
  const map = new Map();
  list.forEach((item) => {
    const key = keyGetter(item);
    const collection = map.get(key);
    if (!collection) {
      map.set(key, [item]);
    } else {
      collection.push(item);
    }
  });
  return map;
}
