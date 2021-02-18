import React, { useMemo } from 'react';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../services/store/publications';
import BTable from 'react-bootstrap/Table';
import { DropdownButton, Form, Row, Col } from 'react-bootstrap';
import { useTable } from 'react-table';
import './publications.scss';

function Table({ title, columns, data, hidden, merge }) {
  const {
    getTableProps,
    headerGroups,
    rows,
    prepareRow,
    setHiddenColumns,
    state: tableState,
  } = useTable({
    columns,
    data,
    initialState: { hiddenColumns: hidden },
  });

  return (
    <div className="mb-4">
      <Row>
        <Col md="8">
          <h3>{title}</h3>
        </Col>
        <Col md="3"></Col>
        <Col md="1">
          <div id={`${title.trim()}-controls`} className="d-flex mb-2">
            <DropdownButton title="Columns">
              <Form>
                {columns.map(({ Header: col }) => {
                  // ignore DOI column
                  if (col != 'DOI')
                    return (
                      <Form.Group
                        controlId={`${col}Visibility`}
                        className="my-1 px-2 "
                      >
                        <Form.Check
                          type="checkbox"
                          label={col}
                          checked={tableState.hiddenColumns.indexOf(col) == -1}
                          onClick={() =>
                            setHiddenColumns((hiddenColumns) => {
                              const index = hiddenColumns.indexOf(col);
                              const newHidden =
                                index > -1
                                  ? hiddenColumns.filter((c) => c != col)
                                  : [...hiddenColumns, col];

                              merge({ hidden: newHidden });
                              return newHidden;
                            })
                          }
                        />
                      </Form.Group>
                    );
                })}
              </Form>
            </DropdownButton>
          </div>
        </Col>
      </Row>

      <BTable responsive striped bordered hover size="sm" {...getTableProps()}>
        <thead>
          {headerGroups.map((headerGroup) => (
            <tr {...headerGroup.getHeaderGroupProps()}>
              {headerGroup.headers.map((column) => (
                <th {...column.getHeaderProps()}>{column.render('Header')}</th>
              ))}
            </tr>
          ))}
        </thead>
        <tbody>
          {rows.map((row, i) => {
            prepareRow(row);
            return (
              <tr {...row.getRowProps()}>
                {row.cells.map((cell) => {
                  return (
                    <td {...cell.getCellProps()}>{cell.render('Cell')}</td>
                  );
                })}
              </tr>
            );
          })}
        </tbody>
      </BTable>
    </div>
  );
}

export default function Publications() {
  // data is loaded in components/app.js
  const state = useSelector((state) => state.publications);
  const tables = useMemo(() => state, [state]) || {};
  const dispatch = useDispatch();

  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <div className="mb-4">
          <h1 className="">Publications</h1>
          <hr />
          <p>
            An overview of published papers, tools, websites and databases
            related to mutational signatures analysis.
          </p>
        </div>

        {tables.orA.data && (
          <Table
            title="Original Research Papers Including Specific Mutational Signatures in mSigPortal"
            columns={tables.orA.columns}
            data={tables.orA.data}
            hidden={tables.orA.hidden}
            merge={(state) => dispatch(actions.mergeState({ orA: state }))}
          />
        )}
        {tables.orB.data && (
          <Table
            title="Original Research Papers involved Mutational Signature Analyses"
            columns={tables.orB.columns}
            data={tables.orB.data}
            hidden={tables.orB.hidden}
            merge={(state) => dispatch(actions.mergeState({ orB: state }))}
          />
        )}
        {tables.rp.data && (
          <Table
            title="Review Papers Focued on Mutational Signatures"
            columns={tables.rp.columns}
            data={tables.rp.data}
            hidden={tables.rp.hidden}
            merge={(state) => dispatch(actions.mergeState({ rp: state }))}
          />
        )}
        {tables.cm.data && (
          <Table
            title="Computational Methods, Databases or Websites for Mutational Signature Analyses"
            columns={tables.cm.columns}
            data={tables.cm.data}
            hidden={tables.cm.hidden}
            merge={(state) => dispatch(actions.mergeState({ cm: state }))}
          />
        )}

        <div>
          <h2 className="font-weight-light">Citations</h2>
          <p>TBA</p>
        </div>

        <p>Last update: 20 JAN 2021.</p>
      </div>
    </div>
  );
}
