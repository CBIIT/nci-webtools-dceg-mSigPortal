import { fileURLToPath, pathToFileURL } from "url";
import { createRequire } from "module";
import minimist from "minimist";
import { initializeSchema } from "./services/utils.js";

// determine if this script was launched from the command line
const isMainModule = process.argv[1] === fileURLToPath(import.meta.url);
const require = createRequire(import.meta.url);

if (isMainModule) {
  const config = require("./config.json");
  const args = minimist(process.argv.slice(2));
  const schemaPath = pathToFileURL(args.schema || "./schema.js");

  const { schema } = await import(schemaPath);
  await initializeSchema(config.database, schema);
  console.log("Initialized all tables");
  process.exit(0);
}
