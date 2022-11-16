import { useEffect, useState, useMemo } from 'react';
import Table from '../../../controls/table/table2';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Description from '../../../controls/description/description';
import { useClusteredQuery } from './apiSlice';
import { actions } from '../../../../services/store/visualization';
import { NavHashLink } from 'react-router-hash-link';

export default function ClusteredIdentification() {
  const store = useSelector((state) => state.visualization);
  const { cluster, projectID } = store.userForm;
  // const { projectID } = store.main;

  const dispatch = useDispatch();
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ clustered: state }));
  const [params, setParams] = useState();

  const { data, error, isFetching } = useClusteredQuery(params, {
    skip: !params,
  });

  // query cluster data
  useEffect(() => {
    if (cluster) {
      const params = {
        userId: projectID,
      };

      setParams(params);
    }
  }, [cluster]);

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
    <div className="bg-white border rounded" style={{ minHeight: '500px' }}>
      <Description
        className="p-3"
        less="This analysis identifies the kataegis events from a VCF file input."
        more={
          <>
            <span>
              Kataegis is a localized substitution hypermutation event, often
              characterized by clusters of C&#x3c;T and/or C&#x3c;G mutations,
              commonly at TpCpN trinucleotides (APOBEC mutations). Click{' '}
              <NavHashLink to="/faq#kataegis">here</NavHashLink> for additional
              information about kataegis.
            </span>
            <p className="mt-3">
              To identify kataegis, input the [Minimum Number of Mutations]
              required for kataegis, the [Maximum Distance] between one mutation
              and the next within a given cluster of mutations being considered
              for kataegis, and a [Chromosome] to be highlighted in the rainfall
              plot. By default, all chromosomes will be shown for the kataegis
              identification.
            </p>
          </>
        }
      />
      <hr />

      <LoadingOverlay active={isFetching} />
      {data && (
        <>
          <Plotly
            className="w-100"
            data={data.rainfall.traces}
            layout={data.rainfall.layout}
            config={data.rainfall.config}
          />
          <p className="p-3">
            The table below is a summary of the kataegis identification based on
            the input parameters. The table can be filtered based on any of the
            columns by entering an appropriate value.
          </p>
          <div className="m-3">
            <Table
              className="border"
              data={data.table}
              columns={columns}
              customOptions={{ hideColumns: true }}
              options={{
                initialState: {
                  hiddenColumns: ['id'],
                },
              }}
              striped
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
