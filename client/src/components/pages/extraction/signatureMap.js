import { useSignatureMapQuery } from './apiSlice';
import Plotly from '../../controls/plotly/plot/plot';
import Table from '../../controls/table/table2';

export default function SignatureMap({ state }) {
  const { id, params, manifest } = state;

  const { data, error } = useSignatureMapQuery(
    {
      plotParams: { userId: id },
      tableParams: {
        id,
        context_type: params.args.context_type,
        signatureMap: manifest.signatureMapJson,
      },
    },
    { skip: !id }
  );
  const { plots, table, multi } = data || {};

  console.log(data);
  console.log(error);

  return (
    <div className="bg-white p-3">
      {table && (
        <Table
          data={table.data}
          columns={table.columns}
          customOptions={{ rowSelectRadio: true }}
          className="border"
          striped
        />
      )}
      {/* {multi && <Plotly data={multi.traces} layout={multi.traces} />} */}
      {plots &&
        plots.map((e) => (
          <div className="border rounded mb-3">
            <Plotly {...e} data={e.traces} />
          </div>
        ))}
    </div>
  );
}
