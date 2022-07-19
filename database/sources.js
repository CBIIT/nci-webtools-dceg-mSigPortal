export const sources = [
  {
    sourcePath: "Association/data.csv",
    table: "association",
    description: "association data",
    columns: [
      "study",
      "strategy",
      "cancer",
      "sample",
      "icgcSpecimenId",
      "icgcDonorId",
      "dataSource",
      "dataType",
      "variableName",
      "variableValue",
      "variableValueType",
    ],
  },

  {
    sourcePath: "Exposure/data.csv",
    table: "exposure",
    description: "exposure data",
    columns: [
      "study",
      "strategy",
      "cancer",
      "organ",
      "sample",
      "signatureSetName",
      "signatureName",
      "exposure",
    ],
  },

  {
    sourcePath: "Seqmatrix/data.csv",
    table: "seqmatrix",
    description: "seqmatrix data",
    columns: [
      "study",
      "cancer",
      "sample",
      "strategy",
      "profile",
      "matrix",
      "mutationType",
      "mutations",
    ],
  },

  {
    sourcePath: "Signature/data.csv",
    table: "signature",
    description: "signature data",
    columns: [
      "source",
      "profile",
      "matrix",
      "signatureSetName",
      "strategy",
      "strandInfo",
      "strand",
      "signatureName",
      "mutationType",
      "contribution",
    ],
  },

  {
    description: "Refresh materialized views",
    type: "postImport",
    callback: async (connection) => {
      await connection.query('REFRESH MATERIALIZED VIEW "seqmatrixOption"');
    }
  }
];
