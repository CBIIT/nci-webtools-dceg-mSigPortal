import { Row, Col, Button, Container } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

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
        ? e.referenceSignature == referenceSignature
        : true)
  )[0];

  // split description string delimited by key:
  const descriptionRegex = /(\w*\s?\w+:)/g;
  const description = metadata?.description
    .split(descriptionRegex)
    .filter((e) => e.length);

  function getStudy() {
    if (data[category] && data[category].length) {
      return [
        ...new Set(
          data[category]
            .filter(
              ({ Etiology, 'Signature Name': signatureName }) =>
                Etiology == etiology && signatureName == signature
            )
            .map((obj) => obj.Study)
        ),
      ].map((Study) => (
        <Col key={Study} lg="2" md="3" sm="4" className="mb-3 d-flex">
          <Button
            size="sm"
            variant="primary"
            onClick={() => mergeEtiology({ study: Study })}
            className={study != Study ? 'disabled' : ''}
            block
          >
            {Study}
          </Button>
        </Col>
      ));
    } else {
      return false;
    }
  }

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
          {metadata.tissueDistribution && (
            <div>
              <strong>Tissue Distribution: </strong>
              {metadata.tissueDistribution}
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
          {description.length && description.length > 1 ? (
            description.map(
              (e, i, arr) =>
                i % 2 == 0 && (
                  <p>
                    <strong>{e}</strong>
                    {arr[i + 1]}
                  </p>
                )
            )
          ) : (
            <p>{description}</p>
          )}

          {/* {profileURL ? (
            <div>
              <SvgContainer
                className="p-3 border rounded mb-3"
                height={'500px'}
                plotPath={profileURL}
                cacheBreaker={false}
              />

              {metadata.Description_strandbias && (
                <p>{metadata.Description_strandbias}</p>
              )}
              {metadata.Description_strandbias && strandbiasURL && (
                <SvgContainer
                  className="p-3 border rounded mb-3"
                  height={'500px'}
                  plotPath={strandbiasURL}
                  cacheBreaker={false}
                />
              )}

              {category == 'Cosmic' && metadata.study && (
                <>
                  <p>
                    Select a cancer study to review the TMB of selected
                    signatures. TMB shown as the numbers of mutations per
                    megabase (log10) attributed to each mutational signature in
                    samples where the signature is present. Only those cancer
                    types with tumors in which signature activity is attributed
                    are shown. The numbers below the dots for each cancer type
                    indicate the number of tumors in which the signatures was
                    attributed (above the horizontal bar, in blue) and the total
                    number of tumors analyzed (below the blue bar, in green).
                  </p>
                  <Row className="justify-content-center">{getStudy()}</Row>
                  {exposureURL.length > 0 && (
                    <SvgContainer
                      className="p-3 border"
                      height={'600px'}
                      plotPath={exposureURL}
                    />
                  )}
                </>
              )}
            </div>
          ) : (
            <div style={{ minHeight: '100px' }}>
              <LoadingOverlay active={true} />
            </div>
          )} */}
        </div>
      )}
    </Container>
  );
}
