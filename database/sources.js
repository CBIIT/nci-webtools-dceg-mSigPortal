import { patternExtractionFormatter } from "./services/formatters.js";

export const sources = [
  {
    sourcePath: "Association/data.csv",
    table: "association",
    description: "association data",
    columns: [
      { sourceName: "Study", name: "study" },
      { sourceName: "Dataset", name: "strategy" },
      { sourceName: "Cancer_Type", name: "cancer" },
      { sourceName: "Sample", name: "sample" },
      { sourceName: "icgc_specimen_id", name: "icgcSpecimenId" },
      { sourceName: "icgc_donor_id", name: "icgcDonorId" },
      { sourceName: "data_source", name: "dataSource" },
      { sourceName: "data_type", name: "dataType" },
      { sourceName: "variable_name", name: "variableName" },
      { sourceName: "variable_value", name: "variableValue" },
      { sourceName: "variable_value_type", name: "variableValueType" },
    ],
  },

  {
    sourcePath: "Exposure/data.csv",
    table: "exposure",
    description: "exposure data",
    columns: [
      { sourceName: "Study", name: "study" },
      { sourceName: "Dataset", name: "strategy" },
      { sourceName: "Cancer_Type", name: "cancer" },
      { sourceName: "Organ", name: "organ" },
      { sourceName: "Sample", name: "sample" },
      { sourceName: "Signature_set_name", name: "signatureSetName" },
      { sourceName: "Signature_name", name: "signatureName" },
      { sourceName: "Exposure", name: "exposure" },
    ],
  },

  {
    sourcePath: "Seqmatrix/data.csv",
    table: "seqmatrix",
    description: "seqmatrix data",
    columns: [
      { sourceName: "Study", name: "study" },
      { sourceName: "Dataset", name: "strategy" },
      { sourceName: "Cancer_Type", name: "cancer" },
      { sourceName: "Sample", name: "sample" },
      { sourceName: "Profile", name: "profile", formatter: patternExtractionFormatter(/^[A-Z]/i) },
      { sourceName: "Profile", name: "matrix", formatter: patternExtractionFormatter(/[0-9]$/) },
      { sourceName: "MutationType", name: "mutationType" },
      { sourceName: "Mutations", name: "mutations" },
    ],
  },
  
  {
    sourcePath: "Signature/data.csv",
    table: "signature",
    description: "signature data",
    columns: [
      { sourceName: "Source", name: "source" },
      { sourceName: "Dataset", name: "strategy" },
      { sourceName: "Profile", name: "profile", formatter: patternExtractionFormatter(/^[A-Z]/i) },
      { sourceName: "Profile", name: "matrix", formatter: patternExtractionFormatter(/[0-9]$/) },
      { sourceName: "Signature_set_name", name: "signatureSetName" },
      { sourceName: "Strand_info", name: "strandInfo" },
      { sourceName: "Strand", name: "strand" },
      { sourceName: "Signature_name", name: "signatureName" },
      { sourceName: "MutationType", name: "mutationType" },
      { sourceName: "Contribution", name: "contribution" },
    ],
  },
];

export default sources;
