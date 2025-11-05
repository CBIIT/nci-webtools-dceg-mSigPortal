import { useEffect, useState } from 'react';
import { Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/catalog';
import {
  useEtiologySignatureQuery,
} from '../etiology/apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../controls/plotly/plot/plot';
import SATSSection from '../etiology/SATSSection';

export default function STSSignatureInfo({ data }) {
  const dispatch = useDispatch();
  const mergeSTS = (state) =>
    dispatch(actions.mergeCatalog({ sts: state }));

  const { category, etiology, signature } =
    useSelector((state) => state.catalog.sts);

  // filter by selected etiology and signature
  const metadata = data.filter(
    (e) =>
      e.etiology == etiology &&
      e.signature == signature
  )[0]?.json;

  // fetch signature plot data
  const [signatureParams, setSignatureParams] = useState(false);
  const {
    data: signaturePlot,
    isFetching: fetchingProfile,
    error: sigPlotError,
  } = useEtiologySignatureQuery(signatureParams, {
    skip: !signatureParams,
  });

  // get profile plot
  useEffect(() => {
    if (metadata) {
      const params = getSignatureParams(metadata.signature);
      if (params) {
        setSignatureParams({
          signatureName: metadata.signature,
          ...params,
        });
      }
    }
  }, [metadata]);

  // Get signature parameters for STS category
  function getSignatureParams(signature) {
    return signature.includes('SBS')
      ? {
          signatureSetName: 'SATS_TS_AACR_GENIE_GRCh37_SBS96',
          profile: 'SBS',
          matrix: '96',
        }
      : signature.includes('DBS')
      ? {
          signatureSetName: 'SATS_TS_AACR_GENIE_GRCh37_DBS78', 
          profile: 'DBS',
          matrix: '78',
        }
      : signature.includes('ID')
      ? {
          signatureSetName: 'SATS_TS_AACR_GENIE_GRCh37_ID83',
          profile: 'ID', 
          matrix: '83',
        }
      : false;
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
            {metadata.genomeBuild && (
              <div>
                <strong>Genome Build: </strong>
                {metadata.genomeBuild}
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
            {signature && (
              <SATSSection selectedSignature={signature} />
            )}
          </div>
        </>
      )}
    </div>
  );
}
