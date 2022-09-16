import { useEffect, useState } from 'react';
import { Row, Col, Button, Container } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';
import { useEtiologyDataQuery, useEtiologyProfileQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../controls/plotly/plot/plot';

const actions = { ...catalogActions, ...modalActions };

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

  const [params, setParams] = useState(false);

  const { data: etiologyData, isFetching } = useEtiologyDataQuery(params, {
    skip: !params,
  });

  useEffect(() => {
    if (metadata) {
      setParams({
        signatureName: metadata.signature,
        study: study.split('_')[0],
        strategy: study.split('_')[1],
      });
    }
  }, [metadata, study]);

  // split description string delimited by key:
  const descriptionRegex = /(\w*\s?\w+:)/g;
  const description = metadata?.description
    ? metadata.description.split(descriptionRegex).filter((e) => e.length)
    : '';

  return (
    <Container fluid="xl" className="p-3">
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
          {metadata.tissueDistribution && (
            <div>
              <strong>Tissue Distribution: </strong>
              {metadata.tissueDistribution}
            </div>
          )}
          <div className="my-4">
            <LoadingOverlay active={isFetching} />
            <Row className="justify-content-center">
              {['TCGA_WES', 'PCAWG_WGS', 'ICGC_WGS', 'GEL_WGS'].map((e) => (
                <Col key={e} lg="2" md="3" sm="4" className="mb-3 d-flex">
                  <Button
                    size="sm"
                    variant="primary"
                    onClick={() => mergeEtiology({ study: e })}
                    className={study != e ? 'disabled' : ''}
                    block
                  >
                    {e}
                  </Button>
                </Col>
              ))}
            </Row>
            {etiologyData?.distributionPlot && (
              <Plotly
                data={etiologyData.distributionPlot.traces}
                layout={etiologyData.distributionPlot.layout}
                config={etiologyData.distributionPlot.config}
              />
            )}
          </div>
          {/* {etiologyData?.data && (
            <div>
              <LoadingOverlay active={isFetching} />
              <Row className="justify-content-center">
                {etiologyData.data.map((e) => (
                  <Col
                    key={e.signatureName}
                    lg="2"
                    md="3"
                    sm="4"
                    className="mb-3 d-flex"
                  >
                    <Button
                      size="sm"
                      variant="primary"
                      onClick={() => mergeEtiology({ study: e })}
                      className={study != e ? 'disabled' : ''}
                      block
                    >
                      {e.signatureName}
                    </Button>
                  </Col>
                ))}
              </Row>
            </div>
          )} */}
        </div>
      )}
    </Container>
  );
}
