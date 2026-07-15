import { createTransport } from "nodemailer";
import { renderTemplate } from "./utils.js";
import { getLogger } from "./logger.js";

const logger = getLogger("notifications");

/**
 * Retrieves the SMTP configuration from the environment.
 * @param {any} env
 * @returns {any} SMTP configuration
 */
export function getSmtpConfig(env) {
  const {
    EMAIL_SMTP_HOST,
    EMAIL_SMTP_PORT,
    EMAIL_SMTP_USER,
    EMAIL_SMTP_PASSWORD,
  } = env;

  let config = {
    host: EMAIL_SMTP_HOST,
    port: EMAIL_SMTP_PORT ? Number.parseInt(EMAIL_SMTP_PORT, 10) : 25,
  };

  if (EMAIL_SMTP_USER && EMAIL_SMTP_PASSWORD) {
    config.auth = {
      user: EMAIL_SMTP_USER,
      pass: EMAIL_SMTP_PASSWORD,
    };
  }

  return config;
}

/**
 * Sends an import completion notification to the webadmin.
 * @param {{ status: string, startTime: string, logs?: string, env?: NodeJS.ProcessEnv }} params
 */
export async function sendImportNotification({
  status,
  startTime,
  logs,
  env = process.env,
}) {
  if (!env.EMAIL_ADMIN || !env.EMAIL_SMTP_HOST) {
    logger.warn(
      "Skipping import notification: EMAIL_ADMIN or EMAIL_SMTP_HOST is not set"
    );
    return;
  }

  const succeeded = status === "COMPLETED";
  const tier = (env.APP_TIER || "").toUpperCase();
  const subject = `[${tier}] mSigPortal Data Import ${succeeded ? "Succeeded" : "Failed"}`;

  const transport = createTransport(getSmtpConfig(env));
  return await transport.sendMail({
    from: env.EMAIL_ADMIN,
    to: env.EMAIL_ADMIN,
    subject,
    html: await renderTemplate("templates/import-notification.html", {
      outcome: succeeded ? "succeeded" : "failed",
      startTime,
      logs,
    }),
  });
}
