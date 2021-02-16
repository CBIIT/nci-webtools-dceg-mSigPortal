import React, { useMemo } from 'react';
import { useSelector } from 'react-redux';
import BTable from 'react-bootstrap/Table';
import { useTable } from 'react-table';

function Table({ title, columns, data }) {
  // Use the state and functions returned from useTable to build your UI
  const { getTableProps, headerGroups, rows, prepareRow } = useTable({
    columns,
    data,
  });

  // Render the UI for your table
  return (
    <div className="mb-3">
      <h2 className="font-weight-light">{title}</h2>
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
        {tables.orA && (
          <Table
            title="Original Research A"
            columns={tables.orA.columns}
            data={tables.orA.data}
          />
        )}

        {tables.orB && (
          <Table
            title="Original Research B"
            columns={tables.orB.columns}
            data={tables.orB.data}
          />
        )}

        {tables.rp && (
          <Table
            title="Review Paper"
            columns={tables.rp.columns}
            data={tables.rp.data}
          />
        )}
        {tables.cm && (
          <Table
            title="Computational Methods"
            columns={tables.cm.columns}
            data={tables.cm.data}
          />
        )}
      </div>
    </div>
  );
}
