import React, { useMemo } from 'react';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../services/store/visualization';
import Table from '../../../components/controls/table/table';

export default function KataegisTable() {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);

  const { kataegisData, pagination, hidden } = visualization.kataegis;

  const columns = useMemo(
    () => [
      {
        id: 'sample',
        accessor: 'sample',
        Header: 'Sample',
        sort: false,
        filterValue: 'asdf',
      },
      {
        id: 'chrom',
        accessor: 'chrom',
        Header: 'Chromosome',
      },
      {
        id: 'start',
        accessor: 'start',
        Header: 'Start',
      },
      {
        id: 'end',
        accessor: 'end',
        Header: 'End',
      },
      {
        id: 'chrom.arm',
        // period delimited keys require a function to be accessed as a string instead of an object
        accessor: (a) => a['chrom.arm'],
        Header: 'Chr. Arm',
      },
      {
        id: 'length',
        accessor: 'length',
        Header: 'Length',
      },
      {
        id: 'number.mut',
        accessor: (a) => a['number.mut'],
        Header: 'Number Mut.',
      },
      {
        id: 'weight.C>X',
        accessor: (a) => a['weight.C>X'],
        Header: 'Weight C>X',
      },
      {
        id: 'confidence',
        accessor: 'confidence',
        Header: 'Confidence',
      },
      {
        id: 'annotation',
        accessor: 'annotation',
        Header: 'Annotation',
      },
      {
        id: 'distanceToTSS',
        accessor: 'distanceToTSS',
        Header: 'Distance to TSS',
      },
      {
        id: 'geneName',
        accessor: 'geneName',
        Header: 'Gene Name',
      },
      {
        id: 'geneID',
        accessor: 'geneID',
        Header: 'Gene ID',
      },
    ],
    []
  );

  const data = useMemo(() => kataegisData, []);

  return (
    <div className="mt-3 px-5" style={{ minWidth: '800px' }}>
      <Table
        title="Kataegis Events"
        columns={columns}
        data={data}
        hidden={hidden}
        pagination={pagination}
        mergeState={(state) =>
          dispatch(actions.mergeVisualization({ kataegis: state }))
        }
      />
    </div>
  );
}
