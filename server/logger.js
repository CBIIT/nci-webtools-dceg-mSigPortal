const path = require('path');
const { createLogger, format, transports } = require('winston');
const { folder, level } = require('./config.json').logs;
require('winston-daily-rotate-file');
const { Console, DailyRotateFile } = transports;

module.exports = new createLogger({
  level: level || 'info',
  format: format.combine(
    format.errors({ stack: true }), // <-- use errors format
    format.colorize(),
    format.timestamp(),
    format.prettyPrint(),
    format.json(),
    format.splat(),
    format.label({ label: '[mSigPortal]' }),
    format.timestamp({
      format: 'YYYY-MM-DD HH:mm:ss',
    }),
    format.printf(({ level, message, timestamp, stack }) => {
      if (stack) {
        // print log trace
        return `[${timestamp}] [${level}] ${message} - ${stack}`;
      }
      return `[${timestamp}] [${level}] ${message}`;
    })
  ),
  transports: [
    new Console(),
    new DailyRotateFile({
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
