import { useState, useEffect } from 'react';
import { Row, Col } from 'react-bootstrap';
import Select from '../../controls/select/selectHookForm';
import { useForm } from 'react-hook-form';
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
      contextType: params.args.context_type,
      denovoSigString,
      decompSigString,
    },
    { skip: !denovoSigString || !decompSigString }
  );
  if (error) {
    console.log('error', error);
  }

  const { denovoPlots, refSigPlots } = data || {};

  const refSigOptions = Object.keys(refSigPlots || {})
    .sort((a, b) => {
      const regex = /\((\w+\.\w+)%\)/;
      const aValue = +a.match(regex)[1];
      const bValue = +b.match(regex)[1];
      return bValue - aValue;
    })
    .map((e) => ({
      label: e,
      value: e,
    }));

  const { control, watch, setValue } = useForm({
    defaultValues: { referenceSignature: '' },
  });
  const { referenceSignature } = watch();

  // select initial table row
  useEffect(() => {
    if (table) {
      const row = table.data[0];
      setDenovo(row['De novo extracted']);
      setDecomposed(row['Global NMF Signatures']);
    }
  }, [table]);

  // select initial reference signature
  useEffect(() => {
    if (
      refSigOptions.length &&
      (!referenceSignature ||
        !refSigOptions.some((e) => e.label === referenceSignature.label))
    )
      setValue('referenceSignature', refSigOptions[0]);
  }, [refSigOptions, referenceSignature]);

  // table options
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
      {denovoPlots &&
        denovoPlots.map((e, i) => (
          <div key={i} className="border rounded mb-3">
            <Plotly {...e} data={e.traces} />
          </div>
        ))}
      {decompSigString !== denovoSigString ? (
        <div className="border rounded p-3">
          <Row>
            <Col lg="auto">
              <Select
                name="referenceSignature"
                label="Reference Signature"
                disabled={!refSigOptions.length}
                options={refSigOptions}
                control={control}
              />
            </Col>
          </Row>
        </div>
      ) : (
        <></>
      )}

      {referenceSignature?.value &&
        Object.keys(refSigPlots).length > 0 &&
        refSigPlots[referenceSignature.value] && (
          <div className="border rounded mt-3">
            <Plotly
              data={refSigPlots[referenceSignature.value].traces}
              layout={refSigPlots[referenceSignature.value].layout}
            />
          </div>
        )}
    </div>
  );
}
