import jStat from 'jstat';

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
  var titles = arr[0].split(/\t|\s+/);
  for (var i = 1; i < arr.length - 1; i++) {
    let data = arr[i].split(/\t|\s+/);

    let dataObject;
    if (data.length === 3) {
      if (!isNumber(data[1]) && isNumber(data[2])) {
        dataObject = {
          sample: data[0],
          value1: data[1],
          value2: data[2],
          titles: [titles[1], titles[2]],
        };
      } else if (isNumber(data[1]) && !isNumber(data[2])) {
        dataObject = {
          sample: data[0],
          value1: data[2],
          value2: data[1],
          titles: [titles[2], titles[1]],
        };
      } else {
        dataObject = {
          sample: data[0],
          value1: data[1],
          value2: data[2],
          titles: [titles[1], titles[2]],
        };
      }
    } else {
      dataObject = {
        sample: data[0],
        value1: data[1],
        titles: [titles[1]],
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

export function linearRegression(x, y) {
  var lr = {};
  var n = y.length;
  var sum_x = 0;
  var sum_y = 0;
  var sum_xy = 0;
  var sum_xx = 0;
  var sum_yy = 0;

  for (var i = 0; i < y.length; i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_xy += x[i] * y[i];
    sum_xx += x[i] * x[i];
    sum_yy += y[i] * y[i];
  }

  lr['sl'] = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
  lr['off'] = (sum_y - lr.sl * sum_x) / n;
  lr['r2'] = Math.pow(
    (n * sum_xy - sum_x * sum_y) /
      Math.sqrt((n * sum_xx - sum_x * sum_x) * (n * sum_yy - sum_y * sum_y)),
    2
  );

  return lr;
}

export function getAvg(grades) {
  const total = grades.reduce((acc, c) => acc + c, 0);
  return total / grades.length;
}

export function round(num, decimalPlaces = 0) {
  if (!Number.isFinite(num)) return num;
  const [coefficient, exponent = 0] = String(num).split('e');
  const shifted = Number(`${coefficient}e${Number(exponent) + decimalPlaces}`);
  if (!Number.isFinite(shifted)) return num;
  const [rounded, roundedExponent = 0] = String(Math.round(shifted)).split('e');
  return Number(`${rounded}e${Number(roundedExponent) - decimalPlaces}`);
}

export function pearson_correlation(x, y, confidence = 0.95) {
  let n = x.length;
  let rho = jStat.corrcoeff(x, y);
  let t = rho * Math.sqrt((n - 2) / (1 - rho ** 2));
  let pvalue = jStat.ttest(t, n - 2, 2);
  let t_critical = jStat.studentt.inv(confidence, n - 2);
  let lower_bound = rho - t_critical * Math.sqrt((1 - rho ** 2) / (n - 2));
  let upper_bound = rho + t_critical * Math.sqrt((1 - rho ** 2) / (n - 2));
  let se = Math.sqrt((1 - rho ** 2) / (n - 2));
  let standard_error = se;
  let z_value = t;
  return {
    rho: rho,
    pvalue: pvalue,
    CI: [lower_bound, upper_bound],
    standard_error: standard_error,
    z_value: z_value,
  };
}

export function spearman_correlation(x, y, confidence = 0.95) {
  let n = x.length;
  let rank_x = jStat.rank(x);
  let rank_y = jStat.rank(y);
  let d_squared = 0;
  for (let i = 0; i < n; i++) {
    d_squared += (rank_x[i] - rank_y[i]) ** 2;
  }
  let rho = 1 - (6 * d_squared) / (n * (n ** 2 - 1));
  let z = Math.asin(rho) * Math.sqrt(n - 3);
  let pvalue = jStat.normal.cdf(-Math.abs(z), 0, 1);
  let z_critical = jStat.normal.inv(confidence + (1 - confidence) / 2, 0, 1);
  let lower_bound = Math.sin(Math.asin(rho) - z_critical / Math.sqrt(n - 3));
  let upper_bound = Math.sin(Math.asin(rho) + z_critical / Math.sqrt(n - 3));
  let se = Math.sqrt((1 - rho ** 2) / (n - 3));
  let standard_error = se;
  let z_value = z;
  return {
    rho: rho,
    pvalue: pvalue,
    CI: [lower_bound, upper_bound],
    standard_error: standard_error,
    z_value: z_value,
  };
}

export function calculateSpearman(x, y) {
  validateCorrelationPairs(x, y);
  // Pearson correlation of average ranks handles ties (including zero exposures).
  // The t test and Fisher-z interval on ranks are asymptotic approximations.
  const stats = calculatePearson(jStat.rank(x), jStat.rank(y));
  return {
    n: stats.n,
    rho: stats.pcorr,
    t: stats.statistic,
    pValue: stats.pValue,
    CILower: stats.CILower,
    CIUpper: stats.CIUpper,
    reason: stats.reason,
    testReason: stats.testReason,
    ciReason: stats.ciReason,
  };
}

function validateCorrelationPairs(x, y) {
  if (x.length !== y.length) {
    throw new Error('Correlation requires the same number of paired values.');
  }
  if (!x.every(Number.isFinite) || !y.every(Number.isFinite)) {
    throw new Error('Correlation requires finite numeric values.');
  }
}

export function calculatePearson(x, y) {
  validateCorrelationPairs(x, y);
  const n = x.length;
  const stats = {
    n,
    pcorr: null,
    statistic: null,
    pValue: null,
    CILower: null,
    CIUpper: null,
    ci: [null, null],
    reason: null,
    testReason: null,
    ciReason: null,
  };
  if (n < 2) {
    stats.reason = 'at least two paired samples are required';
    return stats;
  }
  if (
    x.every((value) => value === x[0]) ||
    y.every((value) => value === y[0])
  ) {
    stats.reason = 'one or both signatures have no variation';
    return stats;
  }
  const xMean = jStat.mean(x);
  const yMean = jStat.mean(y);
  let numerator = 0;
  let xDenominator = 0;
  let yDenominator = 0;
  for (let i = 0; i < n; i++) {
    numerator += (x[i] - xMean) * (y[i] - yMean);
    xDenominator += (x[i] - xMean) ** 2;
    yDenominator += (y[i] - yMean) ** 2;
  }
  // Clamp floating-point overshoots at the correlation boundaries.
  const r = Math.max(
    -1,
    Math.min(1, numerator / Math.sqrt(xDenominator * yDenominator))
  );
  stats.pcorr = r;
  if (n > 2) {
    stats.statistic = r * Math.sqrt((n - 2) / (1 - r ** 2));
    // jStat.ttest expects a sample size and subtracts one internally. Use the
    // distribution directly so the correlation test has exactly n - 2 df.
    stats.pValue =
      Math.abs(r) === 1
        ? 0
        : 2 * jStat.studentt.cdf(-Math.abs(stats.statistic), n - 2);
  } else {
    stats.testReason = 'the t test requires at least three paired samples';
  }
  if (n > 3) {
    // Fisher-z 95% interval, matching the Pearson test used by the plot.
    const z = Math.atanh(r);
    const margin = jStat.normal.inv(0.975, 0, 1) / Math.sqrt(n - 3);
    stats.ci = [Math.tanh(z - margin), Math.tanh(z + margin)];
    [stats.CILower, stats.CIUpper] = stats.ci;
  } else {
    stats.ciReason =
      'the confidence interval requires at least four paired samples';
  }
  return stats;
}

export function extractSubstring(str) {
  // Split the string using underscore as the separator
  const arr = str.split('_');

  // Loop through the array to find the substring
  for (let i = 0; i < arr.length; i++) {
    if (arr[i].indexOf('SBS') !== -1) {
      // If "SBS" is found, return it
      return 'SBS';
    } else if (arr[i].indexOf('DBS') !== -1) {
      // If "DBS" is found, return it
      return 'DBS';
    } else if (arr[i].indexOf('ID') !== -1) {
      // If "ID" is found, return it
      return 'ID';
    }
  }

  // If none of the substrings are found, return null
  return null;
}

export function extractLastWord(str) {
  // Split the string using whitespace as the separator
  const arr = str.split('_');

  // Return the last element of the array
  return arr[arr.length - 1];
}

export function arrayContainsTerms(arr, searchTerms) {
  for (let i = 0; i < searchTerms.length; i++) {
    if (arr.some((item) => item.includes(searchTerms[i]))) {
      return true;
    }
  }
  return false;
}
