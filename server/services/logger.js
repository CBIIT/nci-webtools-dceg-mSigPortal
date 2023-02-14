import util from 'util';
import path from 'path';
import winston from 'winston';
import config from '../config.json' assert { type: 'json' };
import 'winston-daily-rotate-file';

const { createLogger, format, transports } = winston;
const { folder, level } = config.logs;

export default new createLogger({
  level: level || 'info',
  format: format.combine(
    format.colorize(),
    format.splat(),
    format.timestamp({
      format: 'YYYY-MM-DD HH:mm:ss',
    }),
    format.printf(
      (info) =>
        `[${info.timestamp}] [${info.level}] ${
          info.stack || util.format(info.message)
        }`
    )
  ),
  transports: [
    new transports.Console(),
    new transports.DailyRotateFile({
      filename: path.resolve(folder, 'application-%DATE%.log'),
      datePattern: 'YYYY-MM-DD-HH',
      zippedArchive: false,
      maxSize: '1024m',
      timestamp: true,
      maxFiles: '1d',
      prepend: true,
    }),
  ],
  exitOnError: false,
});
