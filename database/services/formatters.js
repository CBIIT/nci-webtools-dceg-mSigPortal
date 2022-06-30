import { format } from "util";

export function ageFormatter(value) {
  if (!isNaN(value)) {
    return value;
  } else if (value) {
    const monthMatches = value.match(/(\d+)m/i);
    return monthMatches ? +monthMatches[1] / 12 : null;
  } else {
    return null;
  }
}

export function patternExtractionFormatter(regex, passthroughIfNoMatch = false) {
  return function (value) {
    if (value === null || value === undefined) return null;
    const defaultValue = passthroughIfNoMatch ? value : null;
    const matches = value.match(regex);
    return matches ? matches[1] : defaultValue;
  };
}

export function chromosomeFormatter(value) {
  if (!value || !value.length) {
    return null;
  } else {
    value = String(value).replace(/^chr/, "").toLowerCase();
    value = { x: 23, y: 24 }[value] || value;
    return +value || null;
  }
}

// this formatter should be removed once the import data format has been fixed
export function unparsedColonNumericValueFormatter(value, record, columnName) {
  if (isNaN(value) || [null, undefined, ""].includes(value)) {
    const colonParser = patternExtractionFormatter(/^.*:(.*)/);
    return colonParser(record[columnName]);
  } else {
    return +value;
  }
}

export function invalidNumberFormatter(value, record, columnName, logger) {
  if ([null, undefined, ""].includes(value)) {
    return null;
  } else if (isNaN(value)) {
    logger?.warn(`Invalid number: ${value} in ${columnName} for ${format(record)}`);
    return null;
  } else {
    return +value;
  }
}
