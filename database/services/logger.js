import util from "util";
import { createLogger, format, transports } from "winston";

export const defaultLogTransports = [new transports.Console()];

export function getLogger(name, level = "info", transports = defaultLogTransports) {
  return new createLogger({
    level,
    format: format.combine(
      format.timestamp({ format: "YYYY-MM-DD HH:mm:ss" }),
      format.label({ label: name }),
      format.printf(({ label, timestamp, level, message }) =>
        [[label, process.pid, timestamp, level].map((s) => `[${s}]`).join(" "), util.format(message)].join(" - ")
      )
    ),
    transports,
    exitOnError: false,
  });
}
