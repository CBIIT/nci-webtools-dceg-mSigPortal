import React, { useMemo } from 'react';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../services/store/publications';
import BTable from 'react-bootstrap/Table';
import { DropdownButton, Form, Row, Col } from 'react-bootstrap';
import { useTable, useGlobalFilter, useAsyncDebounce } from 'react-table';
import './publications.scss';

function GlobalFilter({ globalFilter, setGlobalFilter, handleSearch }) {
  const [value, setValue] = React.useState(globalFilter);
  const onChange = useAsyncDebounce((value) => {
    setGlobalFilter(value || '');
    handleSearch(value || '');
  }, 200);

  return (
    <Form.Group className="m-0">
      <Form.Control
        type="text"
        placeholder="Search"
        value={value || ''}
        onChange={(e) => {
          setValue(e.target.value);
          onChange(e.target.value);
        }}
      />
    </Form.Group>
  );
}

function Table({
  title,
  columns,
  data,
  hidden,
  search,
  handleCheck,
  handleSearch,
}) {
  const {
    getTableProps,
    headerGroups,
    rows,
    prepareRow,
    setHiddenColumns,
    state,
    setGlobalFilter,
  } = useTable(
    {
      columns,
      data,
      initialState: { hiddenColumns: hidden, globalFilter: search },
    },
    useGlobalFilter
  );

  return (
    <div className="mb-4">
      <Row className="mb-2">
        <Col md="8">
          <h3 className="mb-0">{title}</h3>
        </Col>
        <Col md="3">
          <GlobalFilter
            globalFilter={state.globalFilter}
            setGlobalFilter={setGlobalFilter}
            handleSearch={handleSearch}
          />
        </Col>
        <Col md="1">
          <DropdownButton
            title="Columns"
            id={`${title.replace(/\s/g, '')}-controls`}
          >
            <Form>
              {columns.map(({ Header: col }) => {
                // ignore DOI column
                if (col != 'DOI')
                  return (
                    <Form.Group
                      key={`${title.replace(/\s/g, '')}-${col}`}
                      controlId={`${title.replace(
                        /\s/g,
                        ''
                      )}-${col}-visibility`}
                      className="my-1 px-2"
                    >
                      <Form.Check
                        type="checkbox"
                        label={col}
                        checked={state.hiddenColumns.indexOf(col) == -1}
                        onChange={() =>
                          setHiddenColumns((hiddenColumns) => {
                            const index = hiddenColumns.indexOf(col);
                            const newHidden =
                              index > -1
                                ? hiddenColumns.filter((c) => c != col)
                                : [...hiddenColumns, col];

                            handleCheck({ hidden: newHidden });
                            return newHidden;
                          })
                        }
                      />
                    </Form.Group>
                  );
              })}
            </Form>
          </DropdownButton>
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
            search={tables.orA.search}
            handleCheck={(state) =>
              dispatch(actions.mergeState({ orA: state }))
            }
            handleSearch={(e) =>
              dispatch(actions.mergeState({ orA: { search: e } }))
            }
          />
        )}
        {tables.orB.data && (
          <Table
            title="Original Research Papers involved Mutational Signature Analyses"
            columns={tables.orB.columns}
            data={tables.orB.data}
            hidden={tables.orB.hidden}
            search={tables.orB.search}
            handleCheck={(state) =>
              dispatch(actions.mergeState({ orB: state }))
            }
            handleSearch={(e) =>
              dispatch(actions.mergeState({ orB: { search: e } }))
            }
          />
        )}
        {tables.rp.data && (
          <Table
            title="Review Papers Focued on Mutational Signatures"
            columns={tables.rp.columns}
            data={tables.rp.data}
            hidden={tables.rp.hidden}
            search={tables.rp.search}
            handleCheck={(state) => dispatch(actions.mergeState({ rp: state }))}
            handleSearch={(e) =>
              dispatch(actions.mergeState({ rp: { search: e } }))
            }
          />
        )}
        {tables.cm.data && (
          <Table
            title="Computational Methods, Databases or Websites for Mutational Signature Analyses"
            columns={tables.cm.columns}
            data={tables.cm.data}
            hidden={tables.cm.hidden}
            search={tables.cm.search}
            handleCheck={(state) => dispatch(actions.mergeState({ cm: state }))}
            handleSearch={(e) =>
              dispatch(actions.mergeState({ cm: { search: e } }))
            }
          />
        )}

        <div className="mb-4">
          <h3>Citations</h3>
        </div>

        <p>Last update: 20 JAN 2021.</p>
      </div>
    </div>
  );
}
