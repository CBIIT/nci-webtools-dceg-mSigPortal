import { Row, Col, Button, Form } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';

const actions = { ...catalogActions, ...modalActions };
const { Check } = Form;

export default function CancerSpecific() {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeCatalog({ etiology: state }));

  const {
    category,
    etiology,
    tissue,
    refSig,
    all,
    data,
    tissueThumbnails,
    refSigThumbnails,
    exposureURL,
    tissueURL,
    refSigURL,
  } = useSelector((state) => state.catalog.etiology);

  function getCancerEtiology() {
    if (data[category] && data[category].length) {
      return (
        <Row className="justify-content-center">
          {[...new Set(data[category].map((obj) => obj.Etiology))]
            .sort()
            .map((Etiology) => (
              <Col key={Etiology} lg="2" md="3" sm="4" className="mb-3 d-flex">
                <Button
                  size="sm"
                  variant="dark"
                  onClick={() =>
                    mergeEtiology({
                      etiology: Etiology,
                      tissue: '',
                      refSig: '',
                    })
                  }
                  className={etiology != Etiology ? 'disabled' : ''}
                  block
                >
                  {Etiology}
                </Button>
              </Col>
            ))}
        </Row>
      );
    } else {
      return (
        <div style={{ minHeight: '100px' }}>
          <LoadingOverlay active={true} />
        </div>
      );
    }
  }

  function getAllTissues() {
    if (tissueThumbnails.length) {
      return tissueThumbnails.map((tissueSig, index) => {
        return (
          <Col key={index} lg="1" md="2" sm="4" className="mb-2 px-1">
            <div
              className={`sigIcon border rounded ${
                etiology != tissueSig.Etiology
                  ? 'inactive'
                  : tissue == tissueSig['Tissue Specific Signature']
                  ? 'active'
                  : ''
              }`}
              title={`${etiology} - ${tissue['Tissue Specific Signature']}`}
              onClick={() =>
                mergeEtiology({
                  etiology: tissueSig.Etiology,
                  tissue: tissueSig['Tissue Specific Signature'],
                })
              }
            >
              <img
                src={tissueSig.thumbnailURL}
                className="w-100"
                // height="110"
                alt=""
              />
              <strong className="sigLabel">
                {tissueSig['Tissue Specific Signature']}
              </strong>
            </div>
          </Col>
        );
      });
    } else {
      return (
        <div style={{ minHeight: '100px' }}>
          <LoadingOverlay active={true} />
        </div>
      );
    }
  }

  function getTissues() {
    if (tissueThumbnails.length) {
      return tissueThumbnails
        .filter(({ Etiology }) => Etiology == etiology)
        .map((tissueSig, index) => {
          return (
            <Col key={index} md="2" sm="4" className="mb-3">
              <div
                className={`sigIcon border rounded ${
                  tissue == tissueSig['Tissue Specific Signature']
                    ? 'active'
                    : ''
                }`}
                title={`${etiology} - ${tissueSig['Tissue Specific Signature']}`}
                onClick={() =>
                  mergeEtiology({
                    tissue: tissueSig['Tissue Specific Signature'],
                    refSig: '',
                  })
                }
              >
                <img
                  src={tissueSig.thumbnailURL}
                  className="w-100"
                  // height="110"
                  alt=""
                />
                <strong className="sigLabel">
                  {tissueSig['Tissue Specific Signature']}
                </strong>
              </div>
            </Col>
          );
        });
    } else {
      return (
        <div style={{ minHeight: '100px' }}>
          <LoadingOverlay active={true} />
        </div>
      );
    }
  }

  function getRefSig() {
    if (refSigThumbnails.length && tissue) {
      return (
        <Row className="justify-content-center">
          {refSigThumbnails
            .filter(
              (v) =>
                v.Etiology == etiology &&
                v['Tissue Specific Signature'] == tissue
            )
            .map((v, index) => {
              return (
                <Col key={index} md="2" sm="4" className="mb-3">
                  <div
                    className={`sigIcon border rounded ${
                      refSig == v['Ref Signature'] ? 'active' : ''
                    }`}
                    title={`${v['Tissue Specific Signature']} - ${v['Ref Signature']}`}
                    onClick={() =>
                      mergeEtiology({ refSig: v['Ref Signature'] })
                    }
                  >
                    <img
                      src={v.thumbnailURL}
                      className="w-100"
                      // height="110"
                      alt=""
                    />
                    <strong className="sigLabel">
                      {`${v['Ref Signature']} (${v['RefSig Proportion']})`}
                    </strong>
                  </div>
                </Col>
              );
            })}
        </Row>
      );
    } else {
      return <p className="text-center text-muted">Select a Signature</p>;
    }
  }

  function getCancerSpecificInfo() {
    if (data[category] && data[category].length && tissue && refSig) {
      let info = data[category].filter(
        (v) =>
          v['Tissue Specific Signature'] == tissue &&
          v['Ref Signature'] == refSig
      );
      if (info.length) {
        info = info[0];
        return (
          <div>
            <div>
              <strong>Tissue Specific Signature: </strong>
              {tissue}
            </div>
            <div>
              <strong>Reference Signature: </strong>
              {refSig}
            </div>
            <div>
              <strong>RefSig Proportion: </strong>
              {info['RefSig Proportion']}
            </div>
            <div>
              <strong>Study: </strong>
              <a href={info.Study_URL} target="_blank" rel="noreferrer">
                {info.Study}
              </a>
            </div>
            <div>
              <strong>Source: </strong>
              <a href={info.Source_URL} target="_blank" rel="noreferrer">
                {info.Source}
              </a>
            </div>
            {typeof info.Description == 'string' ? (
              <p>{info.Description}</p>
            ) : (
              info.Description.map((text, i) => <p key={i}>{text}</p>)
            )}

            <Plot
              className="p-3 border rounded mb-3"
              height={'500px'}
              plotPath={refSigURL}
            />

            <Plot
              className="p-3 border rounded mb-3"
              height={'500px'}
              plotPath={tissueURL}
            />

            <Plot
              className="p-3 border rounded mb-3"
              height={'600px'}
              plotPath={exposureURL}
            />
          </div>
        );
      } else {
        return (
          <p className="d-flex justify-content-center text-muted">
            Error: No data found for {tissue} and {refSig}
          </p>
        );
      }
    }
  }
  return (
    <div>
      <div className="mb-3">
        <h5 className="separator">Tissue Types</h5>
        {getCancerEtiology()}
      </div>
      <div className="mb-3">
        <h5 className="separator">Tissue Specific Signatures</h5>
        <Row className="justify-content-center mb-3">
          <Col sm="auto">
            <p>
              Select a signature to view more info. Choose{' '}
              <b>Selected Etiology</b> to see signatures in the selected
              etiology, or <b>All Signatures</b> to see signatures for every
              etiology.
            </p>
          </Col>
        </Row>
        <Row className="justify-content-center mb-3">
          <Col sm="auto">
            <Check id="selectedEtiology">
              <Check.Input
                type="radio"
                checked={all == false}
                onChange={() => mergeEtiology({ all: false })}
              />
              <Check.Label className="font-weight-normal">
                Selected Etiology
              </Check.Label>
            </Check>
          </Col>
          <Col sm="auto">
            <Check id="allEtiologies">
              <Check.Input
                type="radio"
                checked={all == true}
                onChange={() => mergeEtiology({ all: true })}
              />
              <Check.Label className="font-weight-normal">
                All Signatures
              </Check.Label>
            </Check>
          </Col>
        </Row>
        <Row className={`justify-content-center ${all ? 'd-none' : ''}`}>
          {getTissues()}
        </Row>
        <Row className={`px-2 justify-content-center ${!all ? 'd-none' : ''}`}>
          {getAllTissues()}
        </Row>
      </div>
      <div className="mb-3">
        <h5 className="separator">Reference Signatures</h5>
        {getRefSig()}
        <div className="p-3">{getCancerSpecificInfo()}</div>
      </div>
    </div>
  );
}
