import React, { useEffect, useMemo, useState } from 'react';
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
    <div className="mb-2">
      <h2 className="font-weight-light">{title}</h2>
      <BTable striped bordered hover size="sm" {...getTableProps()}>
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
  const [state, setState] = useState({ data: {} });
  const mergeState = (obj) => setState({ ...state, ...obj });

  // get publication data
  useEffect(() => {
    const getData = async () => {
      const data = await (await fetch(`api/getPublications`)).json();

      const reducer = (acc, column) => [
        ...acc,
        {
          Header: column,
          accessor: column,
        },
      ];

      mergeState({
        orA: {
          columns: [
            ...new Set(
              ...data['Original Research A'].map((row) => Object.keys(row))
            ),
          ].reduce(reducer, []),
          data: data['Original Research A'],
        },
        orB: {
          columns: [
            ...new Set(
              ...data['Orignal Research B'].map((row) => Object.keys(row))
            ),
          ].reduce(reducer, []),

          data: data['Orignal Research B'],
        },
        rp: {
          columns: [
            ...new Set(...data['Review Paper'].map((row) => Object.keys(row))),
          ].reduce(reducer, []),

          data: data['Review Paper'],
        },
        cm: {
          columns: [
            ...new Set(
              ...data['Computational Methods'].map((row) => Object.keys(row))
            ),
          ].reduce(reducer, []),

          data: data['Computational Methods'],
        },
      });
    };

    getData();
  }, []);

  const tables = useMemo(() => state, [state]);
  console.log(tables);

  return (
    <div className="container bg-white bordered p-3">
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
  );
}
