const { ProvidePlugin } = require('webpack');

module.exports = function override(config, env) {
  // disable source maps since they result in oom errors

  config.resolve.fallback = {
    buffer: require.resolve('buffer/'),
    path: require.resolve('path-browserify'),
    fs: false,
    stream: require.resolve('stream-browserify'),
    assert: require.resolve('assert/'),
  };

  config.plugins.push(
    new ProvidePlugin({
      Buffer: ['buffer', 'Buffer'],
    })
  );

  return config;
};
