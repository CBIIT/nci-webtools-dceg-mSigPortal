export const schema = [
  {
    name: "association",
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
      table.index([
        "study",
        "strategy",
        "cancer",
        "sample",
        "profile",
        "matrix",
      ]);
    },
  },

  {
    name: "seqmatrixOption",
    type: "materializedView",
    dependsOn: ["seqmatrix"],
    schema: (view, connection) => {
      const columns = [
        "study",
        "strategy",
        "cancer",
        "sample",
        "profile",
        "matrix",
      ];
      const query = connection("seqmatrix").distinct(columns);
      view.columns(columns);
      return view.as(query);
    },
  },

  {
    name: "seqmatrixSummary",
    type: "materializedView",
    dependsOn: ["seqmatrix"],
    schema: (view, connection) => {
      const columns = [
        "study",
        "strategy",
        "cancer",
        "sample",
        "profile",
        "matrix",
        "logTotalMutations",
        "meanTotalMutations",
      ];
      const totalCountColumns = [
        ...columns.slice(0, -2),
        connection.raw('sum(mutations) as "totalMutations"'),
      ];
      const summaryColumns = [
        ...columns.slice(0, -3),
        connection.raw(`string_agg(cast(matrix as text), '/') as matrix`),
        connection.raw(`log10(sum("totalMutations")) as "logTotalMutations"`),
        connection.raw(`avg("totalMutations") as "meanTotalMutations"`),
      ];
      const totalCountQuery = connection
        .select(totalCountColumns)
        .from("seqmatrix")
        .groupBy(totalCountColumns.slice(0, -1))
        .orderBy(totalCountColumns.slice(0, -1).map((column) => ({ column })));
      const summaryQuery = connection
        .with("totalCounts", totalCountQuery)
        .select(summaryColumns)
        .from("totalCounts")
        .where("totalMutations", ">", 0)
        .groupBy(summaryColumns.slice(0, -3));

      view.columns(columns);
      return view.as(summaryQuery);
    },
  },

  {
    name: "signature",
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
    schema: (table, connection) => {
      table.increments("id");
      table.string("status");
      table.text("log");
      table.timestamp("createdAt").defaultTo(connection.fn.now());
      table.timestamp("updatedAt").defaultTo(connection.fn.now());
    },
  },
];
