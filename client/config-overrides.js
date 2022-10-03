const { ProvidePlugin } = require('webpack');

module.exports = function override(config, env) {
  // disable source maps since they result in oom errors

  config.resolve.fallback = {
    buffer: require.resolve('buffer/'),
    path: require.resolve('path-browserify'),
    fs: false,
  };

  config.plugins.push(
    new ProvidePlugin({
      Buffer: ['buffer', 'Buffer'],
    })
  );

  return config;
};
