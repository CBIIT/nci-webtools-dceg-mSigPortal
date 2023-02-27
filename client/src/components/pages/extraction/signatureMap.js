import { useSignatureMapQuery } from './apiSlice';

export default function SignatureMap({ state }) {
  const { id, params, manifest } = state;
  const { data } = useSignatureMapQuery(
    {
      id,
      context_type: params.args.context_type,
      signatureMap: manifest.signatureMapJson,
    },
    { skip: !id }
  );

  console.log(data);
  return <pre>{JSON.stringify(data)}</pre>;
}
