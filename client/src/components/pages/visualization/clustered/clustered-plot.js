import { useEffect, useState, useMemo } from 'react';
import Table from '../../../controls/table/table2';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useClusteredQuery } from './apiSlice';

export default function ClusteredPlot() {
  const store = useSelector((state) => state.visualization);
  const { cluster, id } = store.userForm;
  const { sample } = store.clustered;

  const [params, setParams] = useState();

  const { data, error, isFetching } = useClusteredQuery(params, {
    skip: !params,
  });

  // query cluster data
  useEffect(() => {
    if (cluster && id && sample) {
      const params = {
        userId: id,
        sample: sample.value,
      };

      setParams(params);
    }
  }, [cluster, sample]);

  const columns = useMemo(
    () => [
      {
        id: 'sample',
        accessor: 'sample',
        Header: 'Sample',
      },
      {
        id: 'geneId',
        accessor: 'geneId',
        Header: 'Gene ID',
        tooltip:
          'gene ID of the gene in which the kataegis has been identified/nearest gene to the kataegis identified',
      },
      {
        id: 'genome',
        accessor: 'genome',
        Header: 'Genome',
      },
      {
        id: 'mutType',
        accessor: 'mutType',
        Header: 'Mutation Type',
      },
      {
        id: 'chr',
        accessor: 'chr',
        Header: 'Chromosome',
        tooltip: 'The chromosome the kataegis occurs on',
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
        id: 'ref',
        accessor: 'ref',
        Header: 'ref',
      },
      {
        id: 'alt',
        accessor: 'alt',
        Header: 'alt',
      },
      {
        id: 'mutClass',
        accessor: 'mutClass',
        Header: 'Mutation Class',
      },
      {
        id: 'IMDplot',
        accessor: 'IMDplot',
        Header: 'IMDplot',
      },
      {
        id: 'group',
        accessor: 'group',
        Header: 'group',
      },
      {
        id: 'IMD',
        accessor: 'IMD',
        Header: 'IMD',
      },
      {
        id: 'VAF/CCF',
        accessor: 'VAF/CCF',
        Header: 'VAF/CCF',
      },
      {
        id: 'subclass',
        accessor: 'subclass',
        Header: 'subclass',
      },
    ],
    []
  );

  return (
    <div>
      <LoadingOverlay active={isFetching} />
      {data && (
        <>
          <Plotly
            className="w-100"
            data={data.traces}
            layout={data.layout}
            config={data.config}
          />
          <hr />
          <p className="p-3">
            The table below is a summary of the kataegis identification based on
            the input parameters. The table can be filtered based on any of the
            columns by entering an appropriate value.
          </p>
          <div className="m-3">
            <Table
              data={data.table}
              columns={columns}
              customOptions={{
                hideColumns: true,
                download: 'clustered_mutations_identification',
              }}
              options={{
                initialState: {
                  hiddenColumns: ['id'],
                },
              }}
              striped
              bordered
            />
          </div>
        </>
      )}
      {error && (
        <p className="text-center">
          An error has occured. Please check your inputs and try again.
        </p>
      )}
    </div>
  );
}
