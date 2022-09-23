import { useEffect, useState } from 'react';
import { Row, Col, Button, Container } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/catalog';
import {
  useEtiologyDistribtuionQuery,
  useEtiologySignatureQuery,
} from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../controls/plotly/plot/plot';

export default function SignatureInfo({ data }) {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeCatalog({ etiology: state }));

  const { category, etiology, signature, referenceSignature, study } =
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

  const [distributionParams, setDistributionParams] = useState(false);
  const { data: distributionPlot, isFetching: fetchingDistribution } =
    useEtiologyDistribtuionQuery(distributionParams, {
      skip: !distributionParams,
    });

  const [signatureParams, setSignatureParams] = useState(false);
  const { data: profilePlot, isFetching: fetchingProfile } =
    useEtiologySignatureQuery(signatureParams, {
      skip: !signatureParams,
    });

  const studyOptions = [
    {
      label: 'PCAWG WGS',
      value: { study: 'PCAWG', strategy: 'WGS' },
    },
    {
      label: 'TCGA WES',
      value: { study: 'TCGA', strategy: 'WES' },
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
      mergeEtiology({ study: studyOptions[0] });
    }
  }, [metadata]);

  // get distribution plot
  useEffect(() => {
    if (metadata) {
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
          signatureSetName: params.signatureSetName,
          profile: params.profile,
        });
      }
    }
  }, [metadata]);

  // return a valid signatureSetName and profile depending on the selected category and signature
  function getSignatureParams(category, signature) {
    if (category == 'Cosmic') {
      return signature.includes('SBS')
        ? {
            signatureSetName: 'COSMIC_v3.3_Signatures_GRCh37_SBS96',
            profile: 'SBS96',
          }
        : signature.includes('DBS')
        ? {
            signatureSetName: 'COSMIC_v3.3_Signatures_GRCh37_DBS78',
            profile: 'DBS78',
          }
        : signature.includes('ID')
        ? {
            signatureSetName: 'COSMIC_v3.3_Signatures_GRCh37_ID83',
            profile: 'ID83',
          }
        : signature.includes('CN')
        ? {
            signatureSetName: 'COSMIC_v3.3_Signatures_GRCh37_CN48',
            profile: 'CN48',
          }
        : false;
    } else if (category == 'CancerSpecificSignatures_2022') {
      return signature.includes('SBS')
        ? {
            signatureSetName: 'Cancer_Reference_Signatures_2022_GRCh37_SBS96',
            profile: 'SBS96',
          }
        : signature.includes('DBS')
        ? {
            signatureSetName: 'Cancer_Reference_Signatures_2022_GRCh37_DBS78',
            profile: 'DBS78',
          }
        : false;
    } else if (category == 'EnviromentalMutagenesis') {
      return signature.includes('SBS')
        ? {
            signatureSetName: 'Environmental_Mutagen_Signatures_GRCh37_SBS96',
            profile: 'SBS96',
          }
        : signature.includes('DBS')
        ? {
            signatureSetName: 'Environmental_Mutagen_Signatures_GRCh37_DBS78',
            profile: 'DBS78',
          }
        : false;
    } else if (category == 'GeneEdits') {
      return signature.includes('SBS')
        ? {
            signatureSetName: 'Gene_Edits_Signatures_GRCh37_SBS96',
            profile: 'SBS96',
          }
        : false;
    } else if (category == 'CancerSpecificSignatures') {
      return {
        signatureSetName: 'Cancer_Reference_Signatures_GRCh37_SBS96',
        // signatureSetName: 'Organ-specific_Cancer_Signatures_GRCh37_SBS96',
        profile: 'SBS96',
      };
    } else if (category == 'CancerTherapies') {
      return signature.includes('SBS')
        ? {
            signatureSetName: 'Cancer_Therapies_Signatures_GRCh37_SBS96',
            profile: 'SBS96',
          }
        : signature.includes('DBS')
        ? {
            signatureSetName: 'Cancer_Therapies_Signatures_GRCh37_DBS78',
            profile: 'DBS78',
          }
        : false;
    } else if (category == 'Others') {
      return false;
      //   return signature.includes('SBS')
      //     ? {
      //         source: 'Published_signatures',
      //         signatureSetName: 'Organ-specific_Cancer_Signatures_GRCh37_SBS96',
      //         profile: 'SBS96',
      //       }
      //     : false;
    }
  }

  // split description string delimited by key:
  const descriptionRegex = /(\w*\s?\w+:)/g;
  const description = metadata?.description
    ? metadata.description.split(descriptionRegex).filter((e) => e.length)
    : '';

  return (
    <Container fluid className="p-3">
      {metadata && (
        <div>
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
          <div className="my-3 border rounded">
            <LoadingOverlay active={fetchingProfile} />
            {profilePlot ? (
              <Plotly
                data={profilePlot.traces}
                layout={profilePlot.layout}
                config={profilePlot.config}
              />
            ) : (
              <div className="text-center my-4">No data available</div>
            )}
          </div>

          {metadata.tissueDistribution && (
            <div>
              <strong>Tissue Distribution: </strong>
              {metadata.tissueDistribution}
            </div>
          )}
          <div className="my-3 pt-3 border rounded">
            <LoadingOverlay active={fetchingDistribution} />
            <Row className="justify-content-center">
              {studyOptions.map((e) => (
                <Col key={e.label} lg="2" md="3" sm="4" className="mb-3 d-flex">
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
            {distributionPlot?.traces.length ? (
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
      )}
    </Container>
  );
}
