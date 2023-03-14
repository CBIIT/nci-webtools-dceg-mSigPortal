import Table from '../../../components/controls/table/table2';
import { usePublicationsQuery } from './apiSlice';
import './publications.scss';

export default function Publications() {
  const { data, error } = usePublicationsQuery();

  const customOptions = {
    disableColumnSearch: true,
    globalSearch: true,
    hideColumns: true,
  };

  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <div className="mb-4">
          <p>
            An overview of published original research and review papers, tools,
            websites and databases related to mutational signature analyses.
          </p>
        </div>
        {error && 'An error occurred while retrieving publications.'}
        {data && data['Original Research A'] && (
          <div className="mb-5">
            <Table
              title="Original Research Papers Including Specific Mutational Signatures in mSigPortal"
              columns={data['Original Research A'].columns}
              data={data['Original Research A'].data}
              striped
              bordered
              options={{
                initialState: {
                  hiddenColumns: [
                    'diseaseOrPhenotypeOrExposure',
                    'firstAuthor',
                    'lastAuthor',
                    'bioRxivOrPubmedId',
                    'note',
                    'doi',
                  ],
                },
              }}
              customOptions={customOptions}
            />
          </div>
        )}
        {data && data['Original Research B'] && (
          <div className="mb-5">
            <Table
              title="Original Research Papers Involved in Mutational Signature Analyses"
              columns={data['Original Research B'].columns}
              data={data['Original Research B'].data}
              striped
              bordered
              options={{
                initialState: {
                  hiddenColumns: [
                    'diseaseOrPhenotypeOrExposure',
                    'firstAuthor',
                    'lastAuthor',
                    'bioRxivOrPubmedId',
                    'doi',
                  ],
                },
              }}
              customOptions={customOptions}
            />
          </div>
        )}
        {data && data['Review Paper'] && (
          <div className="mb-5">
            <Table
              title="Review Papers Focused on Mutational Signatures"
              columns={data['Review Paper'].columns}
              data={data['Review Paper'].data}
              striped
              bordered
              options={{
                initialState: {
                  hiddenColumns: ['bioRxivOrPubmedId', 'doi'],
                },
              }}
              customOptions={customOptions}
            />
          </div>
        )}
        {data && data['Computational Methods'] && (
          <div className="mb-5">
            <Table
              title="Computational Methods, Tools, Databases or Websites for Mutational Signature Analyses"
              columns={data['Computational Methods'].columns}
              data={data['Computational Methods'].data}
              striped
              bordered
              options={{
                initialState: {
                  hiddenColumns: [
                    'programmingLanguage',
                    'firstAuthor',
                    'lastAuthor',
                    'bioRxivOrPubmedId',
                    'doi',
                    'sourceUrl',
                    'note',
                  ],
                },
              }}
              customOptions={customOptions}
            />
          </div>
        )}

        <div className="mb-4">
          <h3>Citations</h3>
          <p>
            If you use the data from this website, please cite our paper:
            <br />
            Tongwu Zhang, Jian Sang, Alyssa Klein…. Maria Teresa Landi.
            Integrative mutational signature portal (mSigPortal) for cancer
            genomic study (manuscript in preparation).
          </p>
        </div>

        <p>Last update: 15 Nov 2021.</p>
      </div>
    </div>
  );
}
