import { inspect } from 'util';
import {
  createLogger as createWinstonLogger,
  format,
  transports,
} from 'winston';
import pick from 'lodash/pick.js';
import isEmpty from 'lodash/isEmpty.js';

export function formatObject(object) {
  if (object instanceof Error) {
    const errorObject = pick(object, [
      'code',
      'message',
      'stack',
      'stdout',
      'stderr',
    ]);
    return formatObject(errorObject);
  } else if (typeof object === 'string' || typeof object === 'number') {
    return String(object);
  } else if (object === null || object === undefined || isEmpty(object)) {
    return '';
  } else {
    return inspect(object, {
      depth: null,
      compact: true,
      breakLength: Infinity,
    });
  }
}

export function createLogger(name, level = 'info') {
  const { APP_TIER, DATADOG_HOST, DATADOG_API_KEY } = process.env;
  const datadogTransportOptions = {
    host: DATADOG_HOST,
    path: `/api/v2/logs?dd-api-key=${DATADOG_API_KEY}&ddsource=nodejs&service=${APP_TIER}-msigportal-batch`,
    ssl: true,
    format: format.json(),
  };

  return new createWinstonLogger({
    level: level,
    format: format.combine(
      format.timestamp({ format: 'YYYY-MM-DD HH:mm:ss' }),
      format.label({ label: name }),
      format.printf(
        (e) =>
          `[${e.label}] [${e.timestamp}] [${e.level}] - ${formatObject(
            e.message
          )}`
      )
    ),
    transports: [
      new transports.Console(),
      new transports.Http(datadogTransportOptions),
    ],
    exitOnError: false,
  });
}
