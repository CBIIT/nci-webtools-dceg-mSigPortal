import { useSignatureMapQuery } from './apiSlice';

export default function SignatureMap({ state }) {
  const { id, params, manifest } = state;

  const { data } = useSignatureMapQuery(
    {
      id,
      context_type: params.args.context_type,
      signatureMapFile: manifest.signatureMapFile,
      decomposedSignatureFile: manifest.decomposedSignatureFile,
    },
    { skip: !id }
  );

  return <pre>{JSON.stringify(data)}</pre>;
}
