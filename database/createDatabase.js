import { fileURLToPath, pathToFileURL } from "url";
import minimist from "minimist";
import { getConfigFromEnv, initializeSchema } from "./services/utils.js";

// determine if this script was launched from the command line
const isMainModule = process.argv[1] === fileURLToPath(import.meta.url);

if (isMainModule) {
  const config = getConfigFromEnv();
  const args = minimist(process.argv.slice(2));
  const schemaPath = pathToFileURL(args.schema || "./schema.js");

  const { schema } = await import(schemaPath);
  await initializeSchema(config.database, schema);
  console.log("Initialized all tables");
  process.exit(0);
}
