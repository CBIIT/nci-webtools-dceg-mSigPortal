import React from 'react';
import { useSelector } from 'react-redux';
import filterFactory, { textFilter } from 'react-bootstrap-table2-filter';
import { dispatchKataegis } from '../../../services/store';
import {
  Table,
  paginationText,
  paginationSizeSelector,
  paginationButton,
  loadingOverlay,
} from '../../controls/table/table';
import paginationFactory from 'react-bootstrap-table2-paginator';

export default function KataegisTable() {
  const { kataegisData, dataField, filters, order, page, size } = useSelector(
    (state) => state.kataegis
  );

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
      text: 'sample',
      headerTitle: (_) => 'Sample',
      ...commonProps('sample'),
      sort: false,
    },
    {
      dataField: 'chrom',
      text: 'Chr.',
      headerTitle: (_) => 'Chromosome',
      ...commonProps('chrom'),
    },
    {
      dataField: 'start',
      text: 'start',
      headerTitle: (_) => 'Start',
      ...commonProps('start'),
    },
    {
      dataField: 'end',
      text: 'end',
      headerTitle: (_) => 'End',
      ...commonProps('end'),
    },
    {
      dataField: 'chrom.arm',
      text: 'chrom.arm',
      headerTitle: (_) => 'chrom.arm',
      ...commonProps('chrom.arm'),
    },
    {
      dataField: 'length',
      text: 'length',
      headerTitle: (_) => 'Length',
      ...commonProps('length'),
    },
    {
      dataField: 'number.mut',
      text: 'number.mut',
      headerTitle: (_) => 'number.mut',
      ...commonProps('number.mut'),
    },
    {
      dataField: 'weight.C>X',
      text: 'weight.C>X',
      headerTitle: (_) => 'weight.C>X',
      ...commonProps('weight.C>X'),
    },
    {
      dataField: 'confidence',
      text: 'confidence',
      headerTitle: (_) => 'Confidence',
      ...commonProps('confidence'),
    },
    {
      dataField: 'annotation',
      text: 'annotation',
      headerTitle: (_) => 'Annotation',
      ...commonProps('annotation'),
    },
    {
      dataField: 'distanceToTSS',
      text: 'distanceToTSS',
      headerTitle: (_) => 'distanceToTSS',
      ...commonProps('distanceToTSS'),
    },
    {
      dataField: 'geneName',
      text: 'geneName',
      headerTitle: (_) => 'geneName',
      ...commonProps('geneName'),
    },
    {
      dataField: 'geneID',
      text: 'geneID',
      headerTitle: (_) => 'geneID',
      ...commonProps('geneID'),
    },
  ];

  function handleSort(field, order) {
    dispatchKataegis({ dataField: field, order: order });
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
      onPageChange: (page, size) =>
        dispatchKataegis({ page: page, size: size }),
      onSizePerPageChange: (size, page) =>
        dispatchKataegis({ page: page, size: size }),
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
