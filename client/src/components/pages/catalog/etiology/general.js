import { Row, Col, Button, Form } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';

const actions = { ...catalogActions, ...modalActions };
const { Check } = Form;

export default function General({ categories }) {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeCatalog({ etiology: state }));

  const {
    category,
    etiology,
    signature,
    study,
    all,
    data,
    thumbnails,
    profileURL,
    exposureURL,
    strandbiasURL,
  } = useSelector((state) => state.catalog.etiology);

  function getEtiologies() {
    if (data[category] && data[category].length) {
      return (
        <>
          <Row className="justify-content-center">
            {[...new Set(data[category].map((obj) => obj.Etiology))]
              .sort()
              .map((Etiology) => (
                <Col
                  key={Etiology}
                  lg={category == 'GeneEdits' ? '1' : '2'}
                  md="3"
                  sm="4"
                  className="mb-3 d-flex"
                >
                  <Button
                    size="sm"
                    variant="dark"
                    onClick={() =>
                      mergeEtiology({
                        etiology: Etiology,
                        signature: '',
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
        </>
      );
    } else {
      return (
        <div style={{ minHeight: '100px' }}>
          <LoadingOverlay active={true} />
        </div>
      );
    }
  }

  function getAllSignatures() {
    if (thumbnails[category] && thumbnails[category].length) {
      return thumbnails[category].map(
        ({ Etiology, signatureName, thumbnailURL, Source_URL }, index) => (
          <Col key={index} lg="1" md="3" sm="4" className="mb-2 px-1">
            <div
              onClick={() =>
                mergeEtiology({
                  etiology: Etiology,
                  signature: signatureName,
                })
              }
              className={`sigIcon border rounded ${
                etiology != Etiology
                  ? 'inactive'
                  : signatureName == signature
                  ? 'active'
                  : ''
              }`}
              title={`${Etiology} - ${signatureName}`}
            >
              <img
                src={thumbnailURL}
                className="w-100"
                // height="70"
                alt=""
              />
              <strong className="sigLabel" style={{ fontSize: '0.8rem' }}>
                {signatureName}
              </strong>
            </div>
          </Col>
        )
      );
    } else {
      return (
        <div style={{ minHeight: '100px' }}>
          <LoadingOverlay active={true} />
        </div>
      );
    }
  }

  function getSignatures() {
    if (thumbnails[category] && thumbnails[category].length) {
      return thumbnails[category]
        .filter(({ Etiology }) => Etiology == etiology)
        .map(({ signatureName, thumbnailURL, Source_URL }, index) => {
          return (
            <Col key={index} md="2" sm="4" className="mb-3">
              <div
                className={`sigIcon border rounded ${
                  signatureName == signature ? 'active' : ''
                }`}
                title={`${etiology} - ${signatureName}`}
                onClick={() =>
                  mergeEtiology({
                    signature: signatureName,
                  })
                }
              >
                <img
                  src={thumbnailURL}
                  className="w-100"
                  // height="110"
                  alt=""
                />
                <strong className="sigLabel">{signatureName}</strong>
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

  function getInfo() {
    if (data[category] && data[category].length && signature) {
      let info = data[category].filter(
        (signature) =>
          signature['Signature Name'] == signature ||
          signature['Signature'] == signature
      );
      if (info.length) {
        info = info[0];

        return (
          <div>
            <div>
              <strong>Signature Name: </strong>
              {info['Signature Name'] || info['Signature']}
            </div>
            {info['Cell Line'] && (
              <div>
                <strong>Cell Line: </strong>
                {info['Cell Line']}
              </div>
            )}
            {info.Mutagen && (
              <div>
                <strong>Mutagen: </strong>
                {info.Mutagen}
              </div>
            )}
            {info.Treatment && (
              <div>
                <strong>Treatment: </strong>
                {info.Treatment}
              </div>
            )}
            {info.Study && (
              <div>
                <strong>Study: </strong>
                <a href={info.Study_URL} target="_blank" rel="noreferrer">
                  {info.Study}
                </a>
              </div>
            )}
            {info.Source && (
              <div>
                <strong>Source: </strong>
                <a href={info.Source_URL} target="_blank" rel="noreferrer">
                  {info.Source}
                </a>
              </div>
            )}
            {info.URL && (
              <div>
                <strong>Source: </strong>
                <a href={info.URL} target="_blank" rel="noreferrer">
                  {info.URL}
                </a>
              </div>
            )}
            {info.Description &&
              (typeof info.Description == 'string' ? (
                <p>{info.Description}</p>
              ) : (
                info.Description.map((text, i) => <p key={i}>{text}</p>)
              ))}

            {profileURL ? (
              <div>
                <SvgContainer
                  className="p-3 border rounded mb-3"
                  height={'500px'}
                  plotPath={profileURL}
                  cacheBreaker={false}
                />

                {info.Description_strandbias && (
                  <p>{info.Description_strandbias}</p>
                )}
                {info.Description_strandbias && strandbiasURL && (
                  <SvgContainer
                    className="p-3 border rounded mb-3"
                    height={'500px'}
                    plotPath={strandbiasURL}
                    cacheBreaker={false}
                  />
                )}

                {category == 'Cosmic Mutational Signatures (v3.2)' &&
                  info.Study && (
                    <>
                      <p>
                        Select a cancer study to review the TMB of selected
                        signatures. TMB shown as the numbers of mutations per
                        megabase (log10) attributed to each mutational signature
                        in samples where the signature is present. Only those
                        cancer types with tumors in which signature activity is
                        attributed are shown. The numbers below the dots for
                        each cancer type indicate the number of tumors in which
                        the signatures was attributed (above the horizontal bar,
                        in blue) and the total number of tumors analyzed (below
                        the blue bar, in green).
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
            )}
          </div>
        );
      } else {
        return (
          <p className="d-flex justify-content-center text-muted">
            Error: No data found for {signature}
          </p>
        );
      }
    }
  }

  return (
    <div>
      <div className="mb-3">
        <h5 className="separator">
          {categories.filter((e) => e.category == category)[0].etiologyTitle}
        </h5>
        <div>{getEtiologies()}</div>
      </div>
      <div>
        <h5 className="separator">Signatures</h5>

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
          {getSignatures()}
        </Row>
        <Row className={`px-2 justify-content-center ${!all ? 'd-none' : ''}`}>
          {getAllSignatures()}
        </Row>
        <div className="p-3">{getInfo()}</div>
      </div>
    </div>
  );
}
