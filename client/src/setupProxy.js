const { createProxyMiddleware } = require('http-proxy-middleware');

module.exports = function (app) {
  app.use(
    createProxyMiddleware('/extraction', {
      target: 'http://localhost:8332',
      changeOrigin: 'true',
    })
  );
  app.use(createProxyMiddleware('/web', { target: 'http://localhost:8330' }));
  app.use(createProxyMiddleware('/api', { target: 'http://localhost:8330' }));
};
