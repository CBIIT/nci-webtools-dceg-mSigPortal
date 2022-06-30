import Transport from "winston-transport";

export class CustomTransport extends Transport {
  constructor(opts = {}) {
    super(opts);
    this.setHandler(opts.handler);
  }

  setHandler(handler) {
    this.handler = handler || (() => {});
  }

  log(info, callback) {
    setImmediate(() => {
      this.emit("logged", info);
    });

    this.handler(info);
    callback();
  }
}
