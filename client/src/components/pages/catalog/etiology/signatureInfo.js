import { useEffect, useState } from 'react';
import { Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/catalog';
import {
  useEtiologyDistribtuionQuery,
  useEtiologySignatureQuery,
  useEtiologyOrganTableQuery,
} from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../controls/plotly/plot/plot';
import Table from '../../../controls/table/table2';

export default function SignatureInfo({ data }) {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeCatalog({ etiology: state }));

  const { category, etiology, signature, referenceSignature, study, cohort } =
    useSelector((state) => state.catalog.etiology);

  // filter by selected etiology and signature
  const metadata = data.filter(
    (e) =>
      e.etiology == etiology &&
      e.signature == signature &&
      (category == 'CancerSpecificSignatures'
        ? e.json.referenceSignature == referenceSignature
        : true)
  )[0]?.json;

  // fetch etiology distribution plot data
  const [distributionParams, setDistributionParams] = useState(false);
  const { data: distributionPlot, isFetching: fetchingDistribution } =
    useEtiologyDistribtuionQuery(distributionParams, {
      skip: !distributionParams,
    });

  // fetch signature plot data
  const [signatureParams, setSignatureParams] = useState(false);
  const {
    data: signaturePlot,
    isFetching: fetchingProfile,
    error: sigPlotError,
  } = useEtiologySignatureQuery(signatureParams, {
    skip: !signatureParams,
  });

  // fetch reference signature plot data
  const {
    data: refSigPlot,
    isFetching: fetchingRefSig,
    error: refSigPlotError,
  } = useEtiologySignatureQuery(
    referenceSignature == 'RefSig'
      ? {
          signatureName: referenceSignature,
          signatureSetName: 'Cancer_Reference_Signatures_GRCh37_RS32',
          profile: 'RS',
          matrix: '32',
        }
      : {
          signatureName: `Ref.Sig.${referenceSignature?.match(/\d+/) || ''}`,
          signatureSetName: 'Cancer_Reference_Signatures_GRCh37_SBS96',
          profile: 'SBS',
          matrix: '96',
        },
    {
      skip: !referenceSignature,
    }
  );

  // fetch organ table data
  const {
    data: organTable,
    isFetching: fetchingOrgan,
    error: organTableError,
  } = useEtiologyOrganTableQuery(
    { signature: metadata?.signature },
    {
      skip: !metadata || category != 'CancerSpecificSignatures_2022',
    }
  );

  const studyOptions = [
    {
      label: 'TCGA WES',
      value: { study: 'TCGA', strategy: 'WES' },
    },
    {
      label: 'PCAWG WGS',
      value: { study: 'PCAWG', strategy: 'WGS' },
    },
    {
      label: 'ICGC WGS',
      value: { study: 'ICGC-Science-2022', strategy: 'WGS' },
    },
    {
      label: 'GEL WGS',
      value: { study: 'GEL-Science-2022', strategy: 'WGS' },
    },
    {
      label: 'Hartwig WGS',
      value: { study: 'Hartwig-Science-2022', strategy: 'WGS' },
    },
  ];

  // set intial study
  useEffect(() => {
    if (!study) {
      mergeEtiology({ study: studyOptions[1] });
    }
  }, [metadata]);

  // get distribution plot
  useEffect(() => {
    if (
      metadata &&
      (category == 'Cosmic' || category == 'CancerSpecificSignatures_2022')
    ) {
      setDistributionParams({
        signatureName: metadata.signature,
        ...study.value,
      });
    }
  }, [metadata, study]);

  // get profile plot
  useEffect(() => {
    if (metadata) {
      const params = getSignatureParams(category, metadata.signature);
      if (params) {
        setSignatureParams({
          signatureName: metadata.signature,
          ...params,
        });
      }
    }
  }, [metadata]);
  // return a valid signatureSetName and profile depending on the selected category and signature
  // refer to Database/Etiology/signature_refsets_etiology.RData for what parameters to use
  function getSignatureParams(category, signature) {
    if (category == 'Cosmic') {
      return signature.includes('SBS')
        ? {
            signatureSetName: 'COSMIC_v3.3_Signatures_GRCh37_SBS96',
            profile: 'SBS',
            matrix: '96',
          }
        : signature.includes('DBS')
        ? {
            signatureSetName: 'COSMIC_v3.3_Signatures_GRCh37_DBS78',
            profile: 'DBS',
            matrix: '78',
          }
        : signature.includes('ID')
        ? {
            signatureSetName: 'COSMIC_v3.3_Signatures_GRCh37_ID83',
            profile: 'ID',
            matrix: '83',
          }
        : signature.includes('CN')
        ? {
            signatureSetName: 'COSMIC_v3.3_Signatures_GRCh37_CN48',
            profile: 'CN',
            matrix: '48',
          }
        : false;
    } else if (category == 'CancerSpecificSignatures_2022') {
      return signature.includes('SBS')
        ? {
            signatureSetName: 'Cancer_Reference_Signatures_2022_GRCh37_SBS96',
            profile: 'SBS',
            matrix: '96',
          }
        : signature.includes('DBS')
        ? {
            signatureSetName: 'Cancer_Reference_Signatures_2022_GRCh37_DBS78',
            profile: 'DBS',
            matrix: '78',
          }
        : false;
    } else if (category == 'EnviromentalMutagenesis') {
      return signature.includes('SBS')
        ? {
            signatureSetName: 'Environmental_Mutagen_Signatures_GRCh37_SBS96',
            signatureName: metadata.signature.replace(/\s+\(SBS\)/, ''),
            profile: 'SBS',
            matrix: '96',
          }
        : signature.includes('DBS')
        ? {
            signatureSetName: 'Environmental_Mutagen_Signatures_GRCh37_DBS78',
            signatureName: metadata.signature.replace(/\s+\(DBS\)/, ''),
            profile: 'DBS',
            matrix: '78',
          }
        : signature.includes('ID')
        ? {
            signatureSetName: 'Environmental_Mutagen_Signatures_GRCh37_ID29',
            signatureName: metadata.signature.replace(/\s+\(ID\)/, ''),
            profile: 'ID',
            matrix: '29',
          }
        : false;
    } else if (category == 'GeneEdits') {
      return {
        signatureSetName: 'Gene_Edits_Signatures_GRCh37_SBS96',
        signatureName: metadata.signature.replace(/\s+\(SBS\)/, ''),
        profile: 'SBS',
        matrix: '96',
      };
    } else if (category == 'CancerSpecificSignatures') {
      return {
        signatureSetName: 'Organ-specific_Cancer_Signatures_GRCh37_SBS96',
        signatureName: metadata.signature + '%',
        profile: 'SBS',
        matrix: '96',
      };
    } else if (category == 'CancerTherapies') {
      return signature.includes('SBS')
        ? {
            signatureSetName: 'Cancer_Therapies_Signatures_GRCh37_SBS96',
            profile: 'SBS',
            matrix: '96',
          }
        : signature.includes('DBS')
        ? {
            signatureSetName: 'Cancer_Therapies_Signatures_GRCh37_DBS78',
            profile: 'DBS',
            matrix: '78',
          }
        : false;
    } else if (category == 'Others') {
      const sbs192Signatures = [
        'SBS_GA_exp_TSB_WES',
        'SBS_pks_transcStrand_WGS',
        'SBS_BaP_exp_TSB_WGS',
      ];
      if (sbs192Signatures.includes(signature)) {
        return {
          signatureSetName: 'Other_Published_Signatures_GRCh37_SBS192',
          profile: 'SBS',
          matrix: '192',
        };
      } else
        return signature.includes('SBS')
          ? {
              signatureSetName: 'Other_Published_Signatures_%',
              profile: 'SBS',
              matrix: '96',
            }
          : signature.includes('DBS')
          ? {
              signatureSetName: 'Other_Published_Signatures_%',
              profile: 'DBS',
              matrix: '78',
            }
          : signature.includes('ID')
          ? {
              signatureSetName: 'Other_Published_Signatures_%',
              profile: 'ID',
              matrix: '83',
            }
          : false;
    }
  }

  // split description string delimited by key:
  const descriptionRegex = /([A-Z]+[a-z]+\s?\w+:)/g;
  const description = metadata?.description
    ? metadata.description.split(descriptionRegex).filter((e) => e.length)
    : '';

  return (
    <div>
      {metadata && (
        <>
          <div className="p-3">
            <div>
              <strong>Etiology: </strong>
              {metadata.etiology}
            </div>
            <div>
              <strong>Signature: </strong>
              {metadata.signature}
            </div>
            {metadata.referenceSignature && (
              <div>
                <strong>Reference Signature : </strong>
                {metadata.referenceSignature}
              </div>
            )}
            {metadata.cellLine && (
              <div>
                <strong>Cell Line: </strong>
                {metadata.cellLine}
              </div>
            )}
            {metadata.mutagen && (
              <div>
                <strong>Mutagen: </strong>
                {metadata.mutagen}
              </div>
            )}
            {metadata.treatment && (
              <div>
                <strong>Treatment: </strong>
                {metadata.treatment}
              </div>
            )}
            {metadata.refSigProportion && (
              <div>
                <strong>Reference Signature Proportion: </strong>
                {metadata.refSigProportion}
              </div>
            )}
            {metadata.genomeBuild && (
              <div>
                <strong>Genome Build: </strong>
                {metadata.genomeBuild}
              </div>
            )}
            {metadata.cosmic_v3_2 && (
              <div>
                <strong>Identified in COSMICv3.2</strong>
              </div>
            )}
            {metadata.refSig_v1 && (
              <div>
                <strong>Identified in RefSigv1</strong>
              </div>
            )}
            {metadata.cohort && (
              <div>
                <strong>Cohort: </strong>
                {metadata.cohort}
              </div>
            )}
            {metadata.signatureExtractionMethod && (
              <div>
                <strong>Signature Extraction Method: </strong>
                {metadata.signatureExtractionMethod}
              </div>
            )}
            {metadata.tumorType && (
              <div>
                <strong>Tumor Type: </strong>
                {metadata.tumorType}
              </div>
            )}
            {metadata.study && (
              <div>
                <strong>Study: </strong>
                <a href={metadata.studyUrl} target="_blank" rel="noreferrer">
                  {metadata.study}
                </a>
              </div>
            )}
            {metadata.signatureSource && (
              <div>
                <strong>Source: </strong>
                <a href={metadata.sourceUrl} target="_blank" rel="noreferrer">
                  {metadata.signatureSource}
                </a>
              </div>
            )}
            {metadata.source && (
              <div>
                <strong>Source: </strong>
                <a href={metadata.sourceUrl} target="_blank" rel="noreferrer">
                  {metadata.source}
                </a>
              </div>
            )}
            {metadata.descriptionStrandBias && (
              <div>
                <strong>Strand Bias: </strong>
                {metadata.descriptionStrandBias}
              </div>
            )}
            {metadata.note && (
              <div>
                <strong>Note: </strong>
                {metadata.note}
              </div>
            )}
            {description.length > 1 ? (
              description.map(
                (e, i, arr) =>
                  i % 2 == 0 && (
                    <p key={i}>
                      <strong>{e}</strong>
                      {arr[i + 1]}
                    </p>
                  )
              )
            ) : (
              <p>{description}</p>
            )}
            {category == 'CancerSpecificSignatures' && (
              <div className="my-3 border rounded">
                <LoadingOverlay active={fetchingProfile} />
                {refSigPlot && !refSigPlotError ? (
                  <Plotly
                    data={refSigPlot.traces}
                    layout={refSigPlot.layout}
                    config={refSigPlot.config}
                  />
                ) : (
                  <div className="text-center my-4">No data available</div>
                )}
              </div>
            )}
            <div className="my-3 border rounded">
              <LoadingOverlay active={fetchingProfile} />
              {signaturePlot && !sigPlotError ? (
                <Plotly
                  data={signaturePlot.traces}
                  layout={signaturePlot.layout}
                  config={signaturePlot.config}
                />
              ) : (
                <div className="text-center my-4">No data available</div>
              )}
            </div>
          </div>
          {(category == 'Cosmic' ||
            category == 'CancerSpecificSignatures_2022') && (
            <>
              <h5 className="separator my-3">Tissue Distribution</h5>
              <div className="p-3">
                <div>{metadata.tissueDistribution}</div>
                <div className="my-3 pt-3 border rounded">
                  <LoadingOverlay active={fetchingDistribution} />
                  <Row className="justify-content-center">
                    {studyOptions.map((e) => (
                      <Col
                        key={e.label}
                        lg="2"
                        md="3"
                        sm="4"
                        className="mb-3 d-flex"
                      >
                        <Button
                          size="sm"
                          variant="primary"
                          onClick={() => mergeEtiology({ study: e })}
                          className={study.label != e.label ? 'disabled' : ''}
                          block
                        >
                          {e.label}
                        </Button>
                      </Col>
                    ))}
                  </Row>
                  {distributionPlot ? (
                    <Plotly
                      data={distributionPlot.traces}
                      layout={distributionPlot.layout}
                      config={distributionPlot.config}
                    />
                  ) : (
                    <div className="text-center my-4">No data available</div>
                  )}
                </div>
              </div>
            </>
          )}
          {organTable && !organTableError && (
            <>
              <h5 className="separator my-3">Organ-Specific Signatures</h5>
              <div className="p-3">
                <div className="my-3">
                  <LoadingOverlay active={fetchingOrgan} />
                  <Row className="justify-content-center">
                    {organTable.cohortOptions.length &&
                      organTable.cohortOptions.map((e) => (
                        <Col
                          key={e}
                          lg="2"
                          md="3"
                          sm="4"
                          className="mb-3 d-flex"
                        >
                          <Button
                            size="sm"
                            variant="primary"
                            onClick={() => mergeEtiology({ cohort: e })}
                            className={cohort != e ? 'disabled' : ''}
                            block
                          >
                            {e}
                          </Button>
                        </Col>
                      ))}
                    <Col lg="2" md="3" sm="4" className="mb-3 d-flex">
                      <Button
                        size="sm"
                        variant="primary"
                        onClick={() => mergeEtiology({ cohort: '' })}
                        className={cohort != '' ? 'disabled' : ''}
                        block
                      >
                        All Cohorts
                      </Button>
                    </Col>
                  </Row>
                  {organTable?.data ? (
                    <Table
                      data={organTable.data.filter((e) =>
                        cohort ? e.cohort == cohort : true
                      )}
                      columns={organTable.columns}
                      options={{
                        initialState: {
                          hiddenColumns: ['id', 'signature'],
                        },
                      }}
                      customOptions={{
                        hideColumns: true,
                        download: 'etiology_data',
                      }}
                      striped
                      bordered
                    />
                  ) : (
                    <div className="text-center my-4">No data available</div>
                  )}
                </div>
              </div>
            </>
          )}
        </>
      )}
    </div>
  );
}
