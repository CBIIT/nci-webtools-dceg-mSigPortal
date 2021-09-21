const path = require('path');
const { createLogger, format, transports } = require('winston');
const { folder, level } = require('./config.json').logs;
require('winston-daily-rotate-file');

module.exports = new createLogger({
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
          info.stack ||
          (typeof info.message === 'string'
            ? info.message
            : JSON.stringify(info.message))
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
