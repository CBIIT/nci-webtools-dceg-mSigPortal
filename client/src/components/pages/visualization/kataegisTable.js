import React from 'react';
import filterFactory, { textFilter } from 'react-bootstrap-table2-filter';
import {
  Table,
  paginationText,
  paginationSizeSelector,
  paginationButton,
} from '../../controls/table/table';
import paginationFactory from 'react-bootstrap-table2-paginator';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../services/store/visualization';

export default function KataegisTable() {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeKataegis = (state) =>
    dispatch(actions.mergeVisualization({ kataegis: state }));

  const { kataegisData, dataField, order, page, size } = visualization.kataegis;
  const commonProps = (placeholder) => ({
    title: true,
    sort: true,
    onSort: handleSort,
    headerClasses: 'overflow-ellipsis',
    classes: 'overflow-ellipsis',
    filter: textFilter({ placeholder }),
  });

  const columns = [
    {
      dataField: 'sample',
      text: 'Sample',
      headerTitle: (_) => 'Sample',
      ...commonProps('Sample'),
      sort: false,
    },
    {
      dataField: 'chrom',
      text: 'Chromosome',
      headerTitle: (_) => 'Chromosome',
      ...commonProps('Chromosome'),
    },
    {
      dataField: 'start',
      text: 'Start',
      headerTitle: (_) => 'Start',
      ...commonProps('Start'),
    },
    {
      dataField: 'end',
      text: 'End',
      headerTitle: (_) => 'End',
      ...commonProps('End'),
    },
    {
      dataField: 'chrom.arm',
      text: 'Chr. Arm',
      headerTitle: (_) => 'Chromosome Arm.',
      ...commonProps('Chr. Arm'),
    },
    {
      dataField: 'length',
      text: 'Length',
      headerTitle: (_) => 'Length',
      ...commonProps('Length'),
    },
    {
      dataField: 'number.mut',
      text: 'Number Mut.',
      headerTitle: (_) => 'Number Mutations',
      ...commonProps('Number Mut.'),
    },
    {
      dataField: 'weight.C>X',
      text: 'Weight C>X',
      headerTitle: (_) => 'Weight C>X',
      ...commonProps('Weight C>X'),
    },
    {
      dataField: 'confidence',
      text: 'Confidence',
      headerTitle: (_) => 'Confidence',
      ...commonProps('Confidence'),
    },
    {
      dataField: 'annotation',
      text: 'Annotation',
      headerTitle: (_) => 'Annotation',
      ...commonProps('Annotation'),
    },
    {
      dataField: 'distanceToTSS',
      text: 'Distance to TSS',
      headerTitle: (_) => 'Distance to TSS',
      ...commonProps('Distance to TSS'),
    },
    {
      dataField: 'geneName',
      text: 'Gene Name',
      headerTitle: (_) => 'Gene Name',
      ...commonProps('Gene Name'),
    },
    {
      dataField: 'geneID',
      text: 'Gene ID',
      headerTitle: (_) => 'Gene ID',
      ...commonProps('Gene ID'),
    },
  ];

  function handleSort(field, order) {
    mergeKataegis({ dataField: field, order: order });
  }

  const tableProps = {
    keyField: 'index',
    loading: false,
    data: kataegisData.map((d, i) => ({ ...d, index: i })),
    columns: columns,
    defaultSorted: [
      {
        dataField: dataField || 'chrom',
        order: order || 'asc',
      },
    ],
    filter: filterFactory(),
    filterPosition: 'top',
    pagination: paginationFactory({
      page: page || 1,
      sizePerPage: size || 10,
      totalSize: kataegisData.length,
      showTotal: kataegisData.length > 0,
      sizePerPageList: [10, 25, 50, 100],
      paginationTotalRenderer: paginationText(
        'Kataegis event',
        'Kataegis events'
      ),
      sizePerPageRenderer: paginationSizeSelector,
      pageButtonRenderer: paginationButton,
      onPageChange: (page, size) => mergeKataegis({ page: page, size: size }),
      onSizePerPageChange: (size, page) =>
        mergeKataegis({ page: page, size: size }),
    }),
  };

  return (
    kataegisData.length > 0 && (
      <div className="mt-3 px-5" style={{ minWidth: '800px' }}>
        <Table wrapperClasses="table-responsive" {...tableProps} />
      </div>
    )
  );
}
