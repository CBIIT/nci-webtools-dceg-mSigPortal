import { useState, useEffect } from 'react';
import {
  useSignatureMapTableQuery,
  useSignatureMapPlotsQuery,
} from './apiSlice';
import Plotly from '../../controls/plotly/plot/plot';
import Table from '../../controls/table/table2';

export default function SignatureMap({ state }) {
  const { id, params, manifest } = state;
  const [denovoSigString, setDenovo] = useState('');
  const [decompSigString, setDecomposed] = useState('');

  const { data: table, error: tableError } = useSignatureMapTableQuery(
    {
      id,
      context_type: params.args.context_type,
      signatureMap: manifest.signatureMapJson,
    },
    { skip: !id }
  );
  const { data, error } = useSignatureMapPlotsQuery(
    {
      userId: id,
      denovoSigString,
      decompSigString,
    },
    { skip: !denovoSigString || !decompSigString }
  );
  const { plots } = data || {};

  // select inital table row
  useEffect(() => {
    if (table) {
      const row = table.data[0];
      setDenovo(row['De novo extracted']);
      setDecomposed(row['Global NMF Signatures']);
    }
  }, [table]);

  const options = {
    initialState: { selectedRowIds: { 0: true } },
    stateReducer: (newState, action) => {
      // Allow only one row to be selected at a time
      if (action.type === 'toggleRowSelected') {
        newState.selectedRowIds = {
          [action.id]: true,
        };
        setDenovo(table.data[action.id]['De novo extracted']);
        setDecomposed(table.data[action.id]['Global NMF Signatures']);
      }
      return newState;
    },
  };

  return (
    <div className="bg-white p-3">
      {table && (
        <Table
          data={table.data}
          columns={table.columns}
          options={options}
          customOptions={{ rowSelectRadio: true }}
          striped
          bordered
        />
      )}
      {plots &&
        plots.map((e, i) => (
          <div key={i} className="border rounded mb-3">
            <Plotly {...e} data={e.traces} />
          </div>
        ))}
    </div>
  );
}
