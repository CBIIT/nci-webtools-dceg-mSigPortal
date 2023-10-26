export default function MsIndividual_Error(rawData, arg, profile) {
  const error = arg.params_activity.signatureSetName
    ? 'Signature SetName: ' +
      arg.params_activity.signatureSetName +
      ' is not supported in MS Individual'
    : profile === 'Data mismatch'
    ? 'Data mismatch, please try again with different files'
    : 'No data found, please try again';
  return { error };
}
