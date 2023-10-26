function dotProduct(v1, v2) {
  let sum = 0;
  for (let i = 0; i < v1.length; i++) {
    sum += v1[i] * v2[i];
  }
  return sum;
}

function magnitude(v) {
  return Math.sqrt(dotProduct(v, v));
}

function cosineSimilarity(v1, v2) {
  return dotProduct(v1, v2) / (magnitude(v1) * magnitude(v2));
}

module.exports = {
  dotProduct,
  magnitude,
  cosineSimilarity,
};
