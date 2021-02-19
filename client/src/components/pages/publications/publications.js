import React, { useMemo } from 'react';
import { useSelector, useDispatch } from 'react-redux';
import BTable from 'react-bootstrap/Table';
import { DropdownButton, Form, Row, Col } from 'react-bootstrap';
import {
  useTable,
  useGlobalFilter,
  useAsyncDebounce,
  useSortBy,
} from 'react-table';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import {
  faSort,
  faSortUp,
  faSortDown,
} from '@fortawesome/free-solid-svg-icons';
import { actions } from '../../../services/store/publications';
import './publications.scss';

function GlobalFilter({ globalFilter, setGlobalFilter, handleSearch, title }) {
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
        aria-label={`${title.replace(/\s/g, '')}-search`}
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
    useGlobalFilter,
    useSortBy
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
            title={title}
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
                <th {...column.getHeaderProps(column.getSortByToggleProps())}>
                  {column.render('Header')}{' '}
                  <span>
                    {column.isSorted ? (
                      column.isSortedDesc ? (
                        <FontAwesomeIcon
                          className="text-primary"
                          icon={faSortDown}
                        />
                      ) : (
                        <FontAwesomeIcon
                          className="text-primary"
                          icon={faSortUp}
                        />
                      )
                    ) : (
                      <FontAwesomeIcon className="text-primary" icon={faSort} />
                    )}
                  </span>
                </th>
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
  // data is retrieved in components/app.js
  const dispatch = useDispatch();
  const { orA, orB, rp, cm } = useSelector((state) => state.publications);
  const { search: oraSearch, hidden: oraHidden, ...oraTable } = orA;
  const { search: orbSearch, hidden: orbHidden, ...orbTable } = orB;
  const { search: rpSearch, hidden: rpHidden, ...rpTable } = rp;
  const { search: cmSearch, hidden: cmHidden, ...cmTable } = cm;
  const oraMemo = useMemo(() => oraTable, [oraTable]) || {};
  const orbMemo = useMemo(() => orbTable, [orbTable]) || {};
  const rpMemo = useMemo(() => rpTable, [rpTable]) || {};
  const cmMemo = useMemo(() => cmTable, [cmTable]) || {};

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

        {oraMemo.data && (
          <Table
            title="Original Research Papers Including Specific Mutational Signatures in mSigPortal"
            columns={oraMemo.columns}
            data={oraMemo.data}
            hidden={oraHidden}
            search={oraSearch}
            handleCheck={(state) =>
              dispatch(actions.mergeState({ orA: state }))
            }
            handleSearch={(e) =>
              dispatch(actions.mergeState({ orA: { search: e } }))
            }
          />
        )}
        {orbMemo.data && (
          <Table
            title="Original Research Papers involved Mutational Signature Analyses"
            columns={orbMemo.columns}
            data={orbMemo.data}
            hidden={orbHidden}
            search={orbSearch}
            handleCheck={(state) =>
              dispatch(actions.mergeState({ orB: state }))
            }
            handleSearch={(e) =>
              dispatch(actions.mergeState({ orB: { search: e } }))
            }
          />
        )}
        {rpMemo.data && (
          <Table
            title="Review Papers Focued on Mutational Signatures"
            columns={rpMemo.columns}
            data={rpMemo.data}
            hidden={rpHidden}
            search={rpSearch}
            handleCheck={(state) => dispatch(actions.mergeState({ rp: state }))}
            handleSearch={(e) =>
              dispatch(actions.mergeState({ rp: { search: e } }))
            }
          />
        )}
        {cmMemo.data && (
          <Table
            title="Computational Methods, Databases or Websites for Mutational Signature Analyses"
            columns={cmMemo.columns}
            data={cmMemo.data}
            hidden={cmHidden}
            search={cmSearch}
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
