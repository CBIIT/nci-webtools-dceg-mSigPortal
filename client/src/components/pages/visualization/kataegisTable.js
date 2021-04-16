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
        disableSortBy: true,
      },
      {
        id: 'chrom',
        accessor: 'chrom',
        Header: 'Chromosome',
        tooltip: 'the chromosome the kataegis occurs on',
      },
      {
        id: 'start',
        accessor: 'start',
        Header: 'Start',
        tooltip: 'the start position of the kataegis',
      },
      {
        id: 'end',
        accessor: 'end',
        Header: 'End',
        tooltip: 'the end position of the kataegis',
      },
      {
        id: 'chrom.arm',
        // period delimited keys require a function to be accessed as a string instead of an object
        accessor: (a) => a['chrom.arm'],
        Header: 'Chr. Arm',
        tooltip:
          'the arm of the chromosome the kataegis is found on (p represents the short arm, and q represents the long arm)',
      },
      {
        id: 'length',
        accessor: 'length',
        Header: 'Length',
        tooltip: 'how many base pairs the kataegis is',
      },
      {
        id: 'number.mut',
        accessor: (a) => a['number.mut'],
        Header: 'Number Mut.',
        tooltip: 'the number of mutations found in the kataegis identified',
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
        tooltip:
          'what region of the genome this region of kataegis is identified in',
      },
      {
        id: 'distanceToTSS',
        accessor: 'distanceToTSS',
        Header: 'Distance to TSS',
        tooltip:
          'the distance from the kataegis to the transcription start site of the nearest gene',
      },
      {
        id: 'geneName',
        accessor: 'geneName',
        Header: 'Gene Name',
        tooltip:
          'gene in which the kataegis has been identified/nearest gene to the kataegis identified',
      },
      {
        id: 'geneID',
        accessor: 'geneID',
        Header: 'Gene ID',
        tooltip:
          'gene ID of the gene in which the kataegis has been identified/nearest gene to the kataegis identified',
      },
    ],
    []
  );

  const data = useMemo(() => kataegisData, []);

  return (
    <div className="mt-3 px-3" style={{ minWidth: '800px' }}>
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
