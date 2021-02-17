import React, { useMemo } from 'react';
import { useSelector } from 'react-redux';
import BTable from 'react-bootstrap/Table';
import { ButtonGroup, Button } from 'react-bootstrap';
import { useTable } from 'react-table';

function Table({ title, columns, data, hidden }) {
  // Use the state and functions returned from useTable to build your UI
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

  // Render the UI for your table
  return (
    <div className="mb-3">
      <h2 className="font-weight-light">{title}</h2>
      <ButtonGroup className="mb-1">
        {hidden.map((col) => (
          <Button
            variant={
              tableState.hiddenColumns.indexOf(col) > -1
                ? 'outline-secondary'
                : 'primary'
            }
            onClick={() =>
              setHiddenColumns((hiddenColumns) => {
                const index = hiddenColumns.indexOf(col);
                return index > -1
                  ? hiddenColumns.filter((c) => c != col)
                  : [...hiddenColumns, col];
              })
            }
          >
            {col}
          </Button>
        ))}
      </ButtonGroup>
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
  const state = useSelector((state) => state.publications);
  const tables = useMemo(() => state, [state]) || {};

  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <div className="mb-4">
          <h1 className="font-weight-light">Publications</h1>
          <hr />
          <p>
            An overview of published papers, tools, websites and databases
            related to mutational signatures analysis.
          </p>
        </div>

        {tables.orA && (
          <Table
            title="Original Research Papers Including Specific Mutational Signatures in mSigPortal"
            columns={tables.orA.columns}
            data={tables.orA.data}
            hidden={[
              'Disease/Phenotype/Exposure',
              'First author',
              'Last Author',
              'Pubmed ID/BioRxiv',
              'DOI',
            ]}
          />
        )}
        {tables.orB && (
          <Table
            title="Original Research Papers involved Mutational Signature Analyses"
            columns={tables.orB.columns}
            data={tables.orB.data}
            hidden={[
              'Disease/Phenotype/Exposure',
              'First author',
              'Last Author',
              'Pubmed ID/BioRxiv',
              'DOI',
            ]}
          />
        )}
        {tables.rp && (
          <Table
            title="Review Papers Focued on Mutational Signatures"
            columns={tables.rp.columns}
            data={tables.rp.data}
            hidden={['First author', 'Last Author', 'Pubmed ID/BioRxiv', 'DOI']}
          />
        )}
        {tables.cm && (
          <Table
            title="Computational Methods, Databases or Websites for Mutational Signature Analyses"
            columns={tables.cm.columns}
            data={tables.cm.data}
            hidden={['Programming']}
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
