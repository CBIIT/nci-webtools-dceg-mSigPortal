import { unlinkSync, readdirSync, writeFileSync, createReadStream } from 'fs';
import path from 'path';
import { stringify } from 'csv-stringify';
import { groupBy } from 'lodash-es';
import mapValues from 'lodash/mapValues.js';
import { randomUUID } from 'crypto';
import knex from 'knex';
import axios from 'axios';
import FormData from 'form-data';
import { mkdirs, writeJson, readJson } from '../../utils.js';
import { getSignatureData, getSeqmatrixData } from '../../query.js';
import { parseCSV } from '../general.js';

export async function exampleProcessor(
  exampleOutputFolderName,
  inputFolderId,
  randomID,
  env
) {
  const params = {
    args: {
      input_type: 'matrix',
      input_data: 'Samples.txt',
      reference_genome: 'GRCh37',
      exome: 'False',
      context_type: 'DBS78',
      minimum_signatures: '1',
      maximum_signatures: '18',
      nmf_replicates: '500',
      resample: 'True',
      seeds: 'random',
      min_nmf_iterations: '10000',
      max_nmf_iterations: '1000000',
      nmf_test_conv: '10000',
      gpu: 'False',
      stability: '0.8',
      min_stability: '0.2',
      combined_stability: '1',
      allow_stability_drop: 'False',
    },
    signatureQuery: {
      signatureSetName: 'COSMIC_v3_Signatures_GRCh37_DBS78',
      profile: 'DBS',
      matrix: '78',
    },
    seqmatrixQuery: {
      study: 'Sherlock-Lung-232',
      cancer: 'LCINS',
      strategy: 'WGS',
      profile: 'DBS',
      matrix: '78',
    },
    id: exampleOutputFolderName,
    email: '',
    jobName: exampleOutputFolderName,
    form: {
      source: 'user',
      study: { label: 'Sherlock-Lung-232', value: 'Sherlock-Lung-232' },
      cancer: { label: 'LCINS', value: 'LCINS' },
      strategy: { label: 'WGS', value: 'WGS' },
      input_type: { label: 'matrix', value: 'matrix' },
      input_data: 'Samples.txt',
      reference_genome: { label: 'GRCh37', value: 'GRCh37' },
      exome: false,
      signatureSetName: {
        label: 'COSMIC_v3_Signatures_GRCh37_DBS78',
        value: 'COSMIC_v3_Signatures_GRCh37_DBS78',
      },
      signatureName: [{ label: 'all', value: 'all' }],
      extractTool: {
        label: 'SigProfilerExtractor',
        value: 'SigProfilerExtractor',
      },
      matrix_normalization: 'gmm',
      nmf_init: 'random',
      precision: 'single',
      email: '',
      jobName: exampleOutputFolderName,
      gpu: false,
      minimum_signatures: '1',
      maximum_signatures: '18',
      nmf_replicates: '500',
      resample: true,
      seeds: 'random',
      min_nmf_iterations: '10000',
      max_nmf_iterations: '1000000',
      nmf_test_conv: '10000',
      nmf_tolerance: '1e-15',
      stability: '0.8',
      min_stability: '0.2',
      combined_stability: '1',
      context_type: { label: 'DBS78', value: 'DBS78' },
    },
  };

  if (inputFolderId.startsWith('Example_')) {
    const parts = inputFolderId.split('_');
    if (parts.length >= 2) {
      const match = parts[1].match(/([A-Za-z]+)(\d+)/);
      if (match) {
        const profileType = match[1];
        const matrixSize = match[2];

        params.args.context_type = profileType + matrixSize;
        params.signatureQuery.signatureSetName = `COSMIC_v3_Signatures_GRCh37_${profileType}${matrixSize}`;
        params.signatureQuery.profile = profileType;
        params.signatureQuery.matrix = matrixSize;
        params.seqmatrixQuery.profile = profileType;
        params.seqmatrixQuery.matrix = matrixSize;
        //params.id = `Example_${profileType}${matrixSize}_SigProfileExtractor_${randomID}`;
        params.form.signatureSetName.value = `COSMIC_v3_Signatures_GRCh37_${profileType}${matrixSize}`;
        params.form.signatureSetName.label = `COSMIC_v3_Signatures_GRCh37_${profileType}${matrixSize}`;
        params.form.context_type.label = `${profileType}${matrixSize}`;
        params.form.context_type.value = `${profileType}${matrixSize}`;
      } else {
        throw new Error(`Invalid example ID: ${exampleOutputFolderName}`);
      }
    } else {
      throw new Error(`Invalid example ID: ${exampleOutputFolderName}`);
    }
  }
  const dbConnection = knex({
    client: 'postgres',
    connection: {
      host: env.POSTGRES_HOST,
      port: env.POSTGRES_PORT,
      user: env.POSTGRES_USER,
      password: env.POSTGRES_PASS,
      database: env.POSTGRES_DB,
    },
  });

  const { args, signatureQuery, seqmatrixQuery, id, email } = params;
  let paths = await getPaths(params, inputFolderId, randomID, env);

  // const submittedTime = new Date(
  //   (await readJson(paths.statusFile)).submittedAt
  // );

  try {
    // logger.info(id);
    // logger.info(params);
    if (!id) throw new Error('Missing id');

    const inputFolder = path.resolve(env.INPUT_FOLDER, inputFolderId);

    const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);

    await mkdirs([paths.inputFolder, paths.outputFolder]);

    //copy example folder into outputFolder
    //await copy(exampleFolderPath, outputFolder);

    await writeJson(paths.paramsFile, params);
    await writeJson(paths.statusFile, {
      ...(await readJson(paths.statusFile)),
      id,
      status: 'IN_PROGRESS',
    });
    await writeJson(
      paths.manifestFile,
      mapValues(paths, (value) => path.parse(value).base)
    );

    // await uploadWorkingDirectory(inputFolder, outputFolder, id, env);

    const connection = dbConnection;
    const limit = false;

    // query signature data
    const signatureColumns = ['signatureName', 'mutationType', 'contribution'];

    const signatureData = await getSignatureData(
      connection,
      signatureQuery,
      signatureColumns,
      limit
    );

    // transform data into format accepted by SigProfiler
    const groupByType = groupBy(signatureData, (e) => e.mutationType);
    const transposeSignature = Object.values(groupByType).map((signatures) =>
      signatures.reduce((obj, e) => {
        return {
          Type: e.mutationType,
          [e.signatureName]: e.contribution,
          ...obj,
        };
      }, {})
    );

    // write data to tsv file
    const signatureFilePath = path.join(outputFolder, 'signature.tsv');
    stringify(
      transposeSignature,
      {
        header: true,
        delimiter: '\t',
      },
      (error, output) => writeFileSync(signatureFilePath, output)
    );

    ////////////////// query seqmatrix data ////////////////////////
    let seqmatrixData;
    let seqmatrixFilePath;
    let seqmatrixFileName;
    let tsvString;
    if (params.form.source === 'public') {
      const seqmatrixColumns = [
        'sample',
        'profile',
        'matrix',
        'mutationType',
        'mutations',
      ];
      seqmatrixData = await getSeqmatrixData(
        connection,
        seqmatrixQuery,
        seqmatrixColumns,
        limit
      );

      // transform data into format accepted by SigProfiler
      // Extract unique mutation types and samples
      const mutationTypes = [
        ...new Set(seqmatrixData.map((d) => d.mutationType)),
      ];
      const samples = [...new Set(seqmatrixData.map((d) => d.sample))];

      // Initialize the result object with mutation types as keys
      const result = mutationTypes.reduce((acc, mt) => {
        acc[mt] = {};
        return acc;
      }, {});

      // Fill in the result object with mutation counts for each sample
      samples.forEach((s) => {
        mutationTypes.forEach((mt) => {
          const count = seqmatrixData.reduce((acc, d) => {
            if (d.sample === s && d.mutationType === mt) {
              return acc + d.mutations;
            } else {
              return acc;
            }
          }, 0);
          result[mt][s] = count;
        });
      });
      // Write the result to a TSV file using csv-stringify
      const headers = ['MutationType', ...samples];
      const rows = Object.entries(result).map(([mt, counts]) => [
        mt,
        ...samples.map((s) => counts[s] || 0),
      ]);
      const tsvData = rows.map((row) => {
        return {
          MutationType: row[0],
          ...Object.fromEntries(
            samples.map((sample, i) => [sample, row[i + 1] || 0])
          ),
        };
      });

      seqmatrixFilePath = path.join(outputFolder, 'seqmatrix.tsv');
      tsvString = await new Promise((resolve, reject) => {
        stringify(
          tsvData,
          { delimiter: '\t', header: true },
          (error, tsvString) => {
            if (error) {
              reject(error);
            } else {
              resolve(tsvString);
            }
          }
        );
      });

      // Write the TSV string to the file using writeFileSync
      writeFileSync(seqmatrixFilePath, tsvString);
      //logger.info("Result written to signature.tsv");

      seqmatrixFileName = path.basename(seqmatrixFilePath);
    } else {
      seqmatrixFilePath = path.resolve(outputFolder, args.context_type, args.input_data);
    }

    // modify and include parameters
    // const transformArgs = {
    //   ...args,
    //   // input_data: path.join(inputFolder, args.input_data),
    //   input_data:
    //     params.form.source === 'public'
    //       ? seqmatrixFilePath
    //       : path.join(inputFolder, args.input_data),
    //   output: path.join(outputFolder),
    //   signature_database: signatureFilePath,
    // };

    // const cliArgs = Object.entries(transformArgs)
    //   .reduce((params, [key, value]) => [...params, `--${key} ${value}`], [])
    //   .join(' ');
    //
    // logger.info(`[${id}] Run extraction`);
    // await execa(
    //   'python3',
    //   ['services/python/mSigPortal-SigProfilerExtractor.py', cliArgs],
    //   { shell: true }
    // )
    //   .pipeStdout(process.stdout)
    //   .pipeStderr(process.stderr);
    //
    // logger.info('Finished Extraction');

    // import signatures data to database
    const decomposedSignatures = await parseCSV(paths.decomposedSignatureFile);
    const denovoSignatures = await parseCSV(paths.denovoSignatureInput);
    function signatureMapping(e) {
      const { MutationType, ...signatures } = e;
      return Object.entries(signatures).map(([signatureName, mutations]) => ({
        signatureName,
        MutationType,
        mutations,
      }));
    }
    // const transformSignatures = [
    //   ...decomposedSignatures.map(signatureMapping).flat(),
    //   ...denovoSignatures.map(signatureMapping).flat(),
    // ];

    // parse signatureMap csv to JSON
    const signatureMap = await parseCSV(paths.signatureMapFile);
    await writeJson(paths.signatureMapJson, signatureMap);

    // run exploration calculation on denovo and decomposed solutions
    let denovoId, decomposedId;
    try {
      //logger.info(`[${id}] Run Denovo Exploration`);

      const denovoFormData = new FormData();
      if (params.form.source === 'public') {
        denovoFormData.append(
          'matrixFile',
          createReadStream(seqmatrixFilePath)
        );
      } else {
        // denovoFormData.append("matrixFile", createReadStream(paths.matrixFile));
        denovoFormData.append(
          'matrixFile',
          createReadStream(seqmatrixFilePath)
        );
      }
      denovoFormData.append(
        'exposureFile',
        createReadStream(paths.denovoExposureInput)
      );
      denovoFormData.append(
        'signatureFile',
        createReadStream(paths.denovoSignatureInput)
      );

      try {
        const denovoUpload = await axios.post(
          `${env.API_BASE_URL}/api/upload/${randomUUID()}`,
          denovoFormData,
          { headers: denovoFormData.getHeaders() }
        );

        const denovoExploration = await axios.post(
          `${env.API_BASE_URL}/api/submitExploration/${denovoUpload.data.id}`,
          {
            matrixFile: path.parse(seqmatrixFilePath).base,

            exposureFile: path.parse(paths.denovoExposureInput).base,
            signatureFile: path.parse(paths.denovoSignatureInput).base,
          }
        );

        denovoId = denovoExploration.data;
      } catch (error) {
        console.log(error);
      }
    } catch (error) {
      //logger.error("Denovo Exploration Error");
      console.log(error);
      throw error.data;
    }

    try {
      //logger.info(`[${id}] Run Decomposed Exploration`);
      const decomposedFormData = new FormData();

      decomposedFormData.append(
        'matrixFile',
        createReadStream(seqmatrixFilePath)
      );

      decomposedFormData.append(
        'exposureFile',
        createReadStream(paths.decomposedExposureInput)
      );

      decomposedFormData.append(
        'signatureFile',
        createReadStream(paths.decomposedSignatureInput)
      );

      const decomposedUpload = await axios.post(
        `${env.API_BASE_URL}/api/upload/${randomUUID()}`,
        decomposedFormData,
        { headers: decomposedFormData.getHeaders() }
      );

      const decomposedExploration = await axios.post(
        `${env.API_BASE_URL}/api/submitExploration/${decomposedUpload.data.id}`,
        {
          matrixFile: path.parse(seqmatrixFilePath).base,
          exposureFile: path.parse(paths.decomposedExposureInput).base,
          signatureFile: path.parse(paths.decomposedSignatureInput).base,
        }
      );

      decomposedId = decomposedExploration.data;
    } catch (error) {
      //logger.error("Decomposed Exploration Error");

      throw error.data;
    }

    // add exploration ids to manifest
    // await writeJson(paths.manifestFile, {
    //   ...mapValues(paths, (value) => path.parse(value).base),
    //   denovoId,
    //   decomposedId,
    // });
    await writeJson(paths.manifestFile, {
      ...mapValues(paths, (value) => path.parse(value).base),
      denovoId,
      decomposedId,
    });

    // write success status
    await writeJson(paths.statusFile, {
      ...(await readJson(paths.statusFile)),
      status: 'COMPLETED',
    });

    let transformSignatures = [];

    // iterate through each object in decomposedSignatures array
    decomposedSignatures.forEach((decomposedSig) => {
      const { MutationsType, ...mutations } = decomposedSig; // extract MutationsType and values
      Object.entries(mutations).forEach(([signatureName, value]) => {
        transformSignatures.push({
          signatureName,
          MutationType: MutationsType,
          mutations: value,
        });
      });
    });

    // iterate through each object in denovoSignatures array
    denovoSignatures.forEach((denovoSig) => {
      const { MutationsType, ...mutations } = denovoSig; // extract MutationsType and values
      Object.entries(mutations).forEach(([signatureName, value]) => {
        transformSignatures.push({
          signatureName,
          MutationType: MutationsType,
          mutations: value,
        });
      });
    });

    const filePath = path.join(outputFolder, 'local.sqlite3');

    // Check if the file exists
    if (readdirSync(outputFolder).includes('local.sqlite3')) {
      // If the file exists, remove it
      try {
        unlinkSync(filePath);
      } catch (err) {
        console.error('Error deleting file:', err);
      }
    } else {
      console.log('File does not exist. Continuing...');
    }

    // Create a new file with the same name
    const localDb = knex({
      client: 'better-sqlite3',
      connection: {
        filename: filePath,
        createIfNotExist: true,
      },
      useNullAsDefault: true,
      pool: { min: 0, max: 100 },
    });

    await importUserSession(localDb, { signature: transformSignatures });

    // logger.debug(
    //   `Execution Time: ${
    //     (new Date().getTime() - submittedTime.getTime()) / 1000
    //   }`
    // );
  } catch (error) {
    console.error(error);
    await writeJson(paths.statusFile, {
      ...(await readJson(paths.statusFile)),
      status: 'FAILED',
      error: { ...error },
      stopped: new Date(),
    });
  } finally {
    return { id: id, status: 'DONE' };
  }
}

async function getPaths(params, exampleId, randomID, env = process.env) {
  const { id, args } = params;
  const inputFolder = path.resolve(env.INPUT_FOLDER, id);
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  const paramsFile = path.resolve(inputFolder, 'params.json');
  const statusFile = path.resolve(outputFolder, 'status.json');
  const manifestFile = path.resolve(outputFolder, 'manifest.json');
  const databaseFile = path.resolve(outputFolder, 'results.db');

  // map files to be used as input for exploration module
  const solutionsFolder = path.resolve(
    outputFolder,
    args.context_type,
    'Suggested_Solution'
  );

  // const denovoFolder = path.resolve(
  //   solutionsFolder,
  //   `${args.context_type}_De-Novo_Solution`
  // );

  const denovoFolder = path.resolve(solutionsFolder, `De_Novo_Solution`);

  const decomposedFolder = path.resolve(
    solutionsFolder,
    `COSMIC_${args.context_type}_Decomposed_Solution`
  );

  // SigProfilerExtraction log
  const extractionLog = path.resolve(outputFolder, 'JOB_METADATA.txt');

  // matrix file - input for extraction and exploration
  const matrixFile =
    params.form.source === 'public'
      ? ''
      : path.resolve(outputFolder, args.input_data);

  // files for denovo exploration input
  // const denovoExposureInput = path.resolve(
  //   denovoFolder,
  //   "Activities",
  //   `${args.context_type}_De-Novo_Activities_refit.txt`
  // );

  const denovoExposureInput = path.resolve(
    denovoFolder,
    'Activities',
    `De_Novo_Activities_refit.txt`
  );

  // const denovoSignatureInput = path.resolve(
  //   denovoFolder,
  //   "Signatures",
  //   `${args.context_type}_De-Novo_Signatures.txt`
  // );

  const denovoSignatureInput = path.resolve(
    denovoFolder,
    'Signatures',
    `De_Novo_Signatures.txt`
  );

  // files for decomposed exploration input
  // const decomposedExposureInput = path.resolve(
  //   decomposedFolder,
  //   "Activities",
  //   `COSMIC_${args.context_type}_Activities.txt`
  // );

  const decomposedExposureInput = path.resolve(
    decomposedFolder,
    'Activities',
    `COSMIC_${args.context_type}_Activities_refit.txt`
  );

  const decomposedSignatureInput = path.resolve(outputFolder, 'signature.tsv');

  // signature map file
  const signatureMapFile = path.resolve(
    decomposedFolder,
    `De_Novo_map_to_COSMIC_${args.context_type}.csv`
  );

  // signature map file
  const signatureMapJson = path.resolve(
    decomposedFolder,
    `De_Novo_map_to_COSMIC_${args.context_type}.json`
  );

  const decomposedSignatureFile = path.resolve(
    decomposedFolder,
    'Signatures',
    `COSMIC_${args.context_type}_Signatures.txt`
  );

  return {
    inputFolder,
    outputFolder,
    paramsFile,
    statusFile,
    manifestFile,
    databaseFile,
    extractionLog,
    matrixFile,
    denovoExposureInput,
    denovoSignatureInput,
    decomposedExposureInput,
    decomposedSignatureInput,
    signatureMapFile,
    signatureMapJson,
    decomposedSignatureFile,
    contextType: args.context_type,
  };
}

const signatureSchema = [
  {
    name: 'signature',
    schema: (table) => {
      table.increments('id');
      table.string('signatureName');
      table.string('mutationType');
      table.integer('mutations');
    },
    index: (table) => {
      table.index(['signatureName']);
    },
  },
  {
    name: 'signatureOption',
    type: 'materializedView',
    dependsOn: ['signature'],
    create: (connection) => {
      const columns = ['signatureName'];
      return connection.raw(
        [
          `CREATE TABLE signatureOption AS`,
          `SELECT DISTINCT ${columns.join(',')}`,
          `FROM signature`,
        ].join(' ')
      );
    },
    index: (table) => {
      table.index(['signatureName']);
    },
  },
];

async function importUserSession(connection, data, schema = signatureSchema) {
  const tables = signatureSchema.filter((e) => !e.type || e.type === 'table');
  const materializedViews = schema.filter((e) => e.type === 'materializedView');
  const indexedTables = schema.filter((s) => typeof s.index === 'function');

  // create tables
  for (const { name, schema } of tables) {
    await connection.schema.createTable(name, (table) =>
      schema(table, connection)
    );
  }

  // import data
  for (const [tableName, tableData] of Object.entries(data))
    await connection.batchInsert(tableName, tableData, 100);
  // create "materialized" style tables

  for (const { create } of materializedViews) {
    await create(connection);
  }

  // index tables
  for (const { name, index } of indexedTables) {
    await connection.schema.table(name, index);
  }
  return true;
}
