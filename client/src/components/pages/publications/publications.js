import React, { useMemo } from 'react';
import { useSelector, useDispatch } from 'react-redux';
import BTable from 'react-bootstrap/Table';
import { DropdownButton, Form, Row, Col, Button } from 'react-bootstrap';
import {
  useTable,
  useGlobalFilter,
  useAsyncDebounce,
  useSortBy,
  usePagination,
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

function PaginationText({
  from,
  to,
  size,
  singular = 'entry',
  plural = 'entries',
}) {
  return size > 0 ? (
    <span className="react-bootstrap-table-pagination-total ml-2 small text-muted">
      Showing&nbsp;
      {1 + to - from < size && `${from} to ${to} of `}
      {size.toLocaleString()}
      {size === 1 ? ` ${singular}` : ` ${plural || singular + 's'}`}
    </span>
  ) : null;
}

function Table({
  title,
  columns,
  data,
  hidden,
  search,
  pagination,
  mergeState,
}) {
  const {
    getTableProps,
    getTableBodyProps,
    headerGroups,
    prepareRow,
    setHiddenColumns,
    setGlobalFilter,
    page,
    canPreviousPage,
    canNextPage,
    pageCount,
    gotoPage,
    nextPage,
    previousPage,
    setPageSize,
    state: { pageIndex, pageSize, globalFilter, hiddenColumns },
  } = useTable(
    {
      columns,
      data,
      initialState: {
        hiddenColumns: hidden,
        globalFilter: search,
        ...pagination,
      },
    },
    useGlobalFilter,
    useSortBy,
    usePagination
  );

  return (
    <div className="mb-4">
      <Row className="mb-2">
        <Col md="8">
          <h3 className="mb-0">{title}</h3>
        </Col>
        <Col md="3">
          <GlobalFilter
            globalFilter={globalFilter}
            setGlobalFilter={setGlobalFilter}
            handleSearch={(query) => mergeState({ search: query })}
            title={title}
          />
        </Col>
        <Col md="1">
          <DropdownButton
            variant="secondary"
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
                        checked={hiddenColumns.indexOf(col) == -1}
                        onChange={() =>
                          setHiddenColumns((hiddenColumns) => {
                            const index = hiddenColumns.indexOf(col);
                            const newHidden =
                              index > -1
                                ? hiddenColumns.filter((c) => c != col)
                                : [...hiddenColumns, col];

                            mergeState({ hidden: newHidden });
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
                      <FontAwesomeIcon className="text-muted" icon={faSort} />
                    )}
                  </span>
                </th>
              ))}
            </tr>
          ))}
        </thead>
        <tbody {...getTableBodyProps}>
          {page.map((row, i) => {
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
      <Row className="pagination">
        <Col>
          <select
            className="form-control-sm"
            value={pageSize}
            onChange={(e) => {
              setPageSize(Number(e.target.value));
              mergeState({ pagination: { pageSize: Number(e.target.value) } });
            }}
            aria-label="Select a pagination size"
          >
            {[10, 25, 50, 100].map((option) => (
              <option key={option} value={option}>
                {option}
              </option>
            ))}
          </select>
          <PaginationText
            from={pageIndex * pageSize + 1}
            to={pageIndex * pageSize + page.length}
            size={data.length}
          />
        </Col>
        <Col className="d-flex">
          <div className="ml-auto">
            {canPreviousPage && (
              <Button
                className="p-2"
                variant="link"
                onClick={() => {
                  gotoPage(0);
                  mergeState({ pagination: { pageIndex: 0 } });
                }}
              >
                {'<<'}
              </Button>
            )}
            {canPreviousPage && (
              <Button
                className="p-2"
                variant="link"
                onClick={() => {
                  previousPage();
                  mergeState({ pagination: { pageIndex: pageIndex - 1 } });
                }}
              >
                {'<'}
              </Button>
            )}
            {canNextPage && (
              <Button
                className="p-2"
                variant="link"
                onClick={() => {
                  nextPage();
                  mergeState({ pagination: { pageIndex: pageIndex + 1 } });
                }}
              >
                {'>'}
              </Button>
            )}
            {canNextPage && (
              <Button
                className="p-2"
                variant="link"
                onClick={() => {
                  gotoPage(pageCount - 1);
                  mergeState({ pagination: { pageIndex: pageCount - 1 } });
                }}
              >
                {'>>'}
              </Button>
            )}
          </div>
        </Col>
      </Row>
    </div>
  );
}

export default function Publications() {
  // data is retrieved in components/app.js
  const dispatch = useDispatch();
  const { orA, orB, rp, cm } = useSelector((state) => state.publications);
  const {
    search: oraSearch,
    hidden: oraHidden,
    pagination: oraPagination,
    ...oraTable
  } = orA;
  const {
    search: orbSearch,
    hidden: orbHidden,
    pagination: orbPagination,
    ...orbTable
  } = orB;
  const {
    search: rpSearch,
    hidden: rpHidden,
    pagination: rpPagination,
    ...rpTable
  } = rp;
  const {
    search: cmSearch,
    hidden: cmHidden,
    pagination: cmPagination,
    ...cmTable
  } = cm;
  const oraMemo = useMemo(() => oraTable, [oraTable]) || {};
  const orbMemo = useMemo(() => orbTable, [orbTable]) || {};
  const rpMemo = useMemo(() => rpTable, [rpTable]) || {};
  const cmMemo = useMemo(() => cmTable, [cmTable]) || {};

  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <div className="mb-4">
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
            pagination={oraPagination}
            mergeState={(state) => dispatch(actions.mergeState({ orA: state }))}
          />
        )}
        {orbMemo.data && (
          <Table
            title="Original Research Papers involved Mutational Signature Analyses"
            columns={orbMemo.columns}
            data={orbMemo.data}
            hidden={orbHidden}
            search={orbSearch}
            pagination={orbPagination}
            mergeState={(state) => dispatch(actions.mergeState({ orB: state }))}
          />
        )}
        {rpMemo.data && (
          <Table
            title="Review Papers Focued on Mutational Signatures"
            columns={rpMemo.columns}
            data={rpMemo.data}
            hidden={rpHidden}
            search={rpSearch}
            pagination={rpPagination}
            mergeState={(state) => dispatch(actions.mergeState({ rp: state }))}
          />
        )}
        {cmMemo.data && (
          <Table
            title="Computational Methods, Databases or Websites for Mutational Signature Analyses"
            columns={cmMemo.columns}
            data={cmMemo.data}
            hidden={cmHidden}
            search={cmSearch}
            pagination={cmPagination}
            mergeState={(state) => dispatch(actions.mergeState({ cm: state }))}
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
