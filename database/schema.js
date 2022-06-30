export const schema = [
  {
    name: "association",
    recreate: true,
    schema: (table) => {
      table.increments("id");
      table.string("study");
      table.string("strategy");
      table.string("cancer");
      table.string("sample");
      table.string("icgcSpecimenId");
      table.string("icgcDonorId");
      table.string("dataSource");
      table.string("dataType");
      table.string("variableName");
      table.string("variableValue");
      table.string("variableValueType");
      table.index(["study", "strategy", "cancer", "sample"]);
    },
  },

  {
    name: "exposure",
    recreate: true,
    schema: (table) => {
      table.increments("id");
      table.string("study");
      table.string("strategy");
      table.string("cancer");
      table.string("organ");
      table.string("sample");
      table.string("signatureSetName");
      table.string("signatureName");
      table.double("exposure");
      table.index(["study", "strategy", "cancer", "signatureSetName"]);
    },
  },

  {
    name: "seqmatrix",
    recreate: true,
    schema: (table) => {
      table.increments("id");
      table.string("study");
      table.string("strategy");
      table.string("cancer");
      table.string("sample");
      table.string("profile");
      table.integer("matrix");
      table.string("mutationType");
      table.integer("mutations");
      table.index(["study", "strategy", "cancer", "sample", "profile", "matrix"]);
    },
  },

  {
    name: "signature",
    recreate: true,
    schema: (table) => {
      table.increments("id");
      table.string("source");
      table.string("strategy");
      table.string("profile");
      table.integer("matrix");
      table.string("signatureSetName");
      table.string("strandInfo");
      table.string("strand");
      table.string("signatureName");
      table.string("mutationType");
      table.double("contribution");
      table.index(["profile", "matrix", "signatureSetName"]);
    },
  },

   {
    name: "importLog",
    recreate: true,
    schema: (table, connection) => {
      table.increments("id");
      table.string("status");
      table.text("log");
      table.boolean("cancelled").defaultTo(false);
      table.timestamp("createdAt").defaultTo(connection.fn.now());
      table.timestamp("updatedAt").defaultTo(connection.fn.now());
    },
  },
];
