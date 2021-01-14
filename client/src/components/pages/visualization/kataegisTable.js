import React from 'react';
import { useSelector } from 'react-redux';
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
  const { kataegisData, dataField, order, page, size } = useSelector(
    (state) => state.kataegis
  );

  const commonProps = {
    title: true,
    sort: true,
    onSort: handleSort,
    // headerStyle: { width: '65px', minWidth: '65px' },
    headerClasses: 'overflow-ellipsis',
    classes: 'overflow-ellipsis',
  };

  const columns = [
    {
      dataField: 'sample',
      text: 'sample',
      headerTitle: (_) => 'Sample',
      ...commonProps,
      sort: false,
    },
    {
      dataField: 'chrom',
      text: 'Chr.',
      headerTitle: (_) => 'Chromosome',
      ...commonProps,
    },
    {
      dataField: 'start',
      text: 'start',
      headerTitle: (_) => 'Start',
      ...commonProps,
    },
    {
      dataField: 'end',
      text: 'end',
      headerTitle: (_) => 'End',
      ...commonProps,
    },
    {
      dataField: 'chrom.arm',
      text: 'chrom.arm',
      headerTitle: (_) => 'chrom.arm',
      ...commonProps,
    },
    {
      dataField: 'length',
      text: 'length',
      headerTitle: (_) => 'Length',
      ...commonProps,
    },
    {
      dataField: 'number.mut',
      text: 'number.mut',
      headerTitle: (_) => 'number.mut',
      ...commonProps,
    },
    {
      dataField: 'weight.C>X',
      text: 'weight.C>X',
      headerTitle: (_) => 'weight.C>X',
      ...commonProps,
    },
    {
      dataField: 'confidence',
      text: 'confidence',
      headerTitle: (_) => 'Confidence',
      ...commonProps,
    },
    {
      dataField: 'annotation',
      text: 'annotation',
      headerTitle: (_) => 'Annotation',
      ...commonProps,
    },
    {
      dataField: 'distanceToTSS',
      text: 'distanceToTSS',
      headerTitle: (_) => 'distanceToTSS',
      ...commonProps,
    },
    {
      dataField: 'geneName',
      text: 'geneName',
      headerTitle: (_) => 'geneName',
      ...commonProps,
    },
    {
      dataField: 'geneID',
      text: 'geneID',
      headerTitle: (_) => 'geneID',
      ...commonProps,
    },
  ];

  function handleSort(field, order) {
    dispatchKataegis({ dataField: field, order: order });
  }

  const tableProps = {
    keyField: 'chrom',
    loading: false,
    data: kataegisData,
    columns: columns,
    overlay: loadingOverlay,
    defaultSorted: [
      {
        dataField: dataField || 'chrom',
        order: order || 'asc',
      },
    ],
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
        <div
          key="controls"
          className="d-flex align-items-center justify-content-between"
        >
          {/* <div key="snpSearch" className="d-flex mb-2">
          <input
            style={{ minWidth: '280px' }}
            className="form-control form-control-sm"
            placeholder="Search for a SNP or SNPs (ex/ rs3 rs4)"
            value={summarySnpTables.snp}
            onChange={(e) => setSnp(e.target.value)}
            aria-label="Filter SNP"
          />
          <button
            className="btn btn-sm btn-silver flex-shrink-auto d-flex"
            onClick={handleSnpReset}
          >
            <Icon className="opacity-50" name="times" width="12" />
            <span className="sr-only">Clear</span>
          </button>
          <button
            className="btn btn-sm btn-silver flex-shrink-auto mx-2"
            onClick={handleSnpLookup}
          >
            Search
          </button>
        </div>
      </div> */}
        </div>
        <Table wrapperClasses="table-responsive" {...tableProps} />
      </div>
    )
  );
}
