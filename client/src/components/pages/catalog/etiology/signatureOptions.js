import { Row, Col, Form } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useThumbnailsQuery } from './apiSlice';

const actions = { ...catalogActions, ...modalActions };
const { Check } = Form;

export default function SignatureOptions({ data }) {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeCatalog({ etiology: state }));

  const {
    showAllSignatures,
    category,
    etiology,
    signature,
    referenceSignature,
  } = useSelector((state) => state.catalog.etiology);

  const {
    data: thumbnails,
    error,
    isFetching,
  } = useThumbnailsQuery(
    {
      keys: data.map((e) => getThumbnailKey(e.signature)),
    },
    { skip: !data }
  );

  // for category CancerSpecificSignatures
  const { data: refSigThumbnails } = useThumbnailsQuery(
    {
      keys: data.map((e) => getThumbnailKey(e.json.referenceSignature)),
    },
    { skip: !data || category != 'CancerSpecificSignatures' }
  );

  // creates s3 path key from signature name
  function getThumbnailKey(signature) {
    // replace slashes with colons in filenames
    const file = signature ? signature.replaceAll(/\//g, ':') : '';
    return encodeURIComponent(
      category == 'CancerSpecificSignatures_2022'
        ? `msigportal/Database/Etiology/Profile_logo/Cancer_Specific_Signature_2022/${file}.svg`
        : `msigportal/Database/Etiology/Profile_logo/${file}.svg`
    );
  }

  function naturalSort(a, b) {
    const sigA = a.signature || false;
    const sigB = b.signature || false;
    if (sigA)
      return sigA.localeCompare(sigB, undefined, {
        numeric: true,
        sensitivity: 'base',
      });
    else {
      return 0;
    }
  }

  function profileSort(a, b) {
    const sigOrder = [/SBS/g, /DBS/g, /ID/g, /CN/g];
    let c = 0,
      d = 0;

    sigOrder.forEach((profile, i) => {
      const sigA = a.signature;
      const sigB = b.signature;
      if (sigA.match(profile)) c = i;
      if (sigB.match(profile)) d = i;
    });

    return c - d;
  }

  // get array with only unique signatures
  const unqiueSignatures = data
    ? data.filter((e, i, arr) => {
        if (i == 0) {
          return e;
        } else if (e.signature !== arr[i - 1].signature) {
          return e;
        }
      })
    : [];

  return (
    <div>
      <h5 className="separator">Signatures</h5>
      <Row className="justify-content-center mb-3">
        <Col sm="auto">
          <p>
            Select a signature to view more info. Choose{' '}
            <b>Selected Etiology</b> to see signatures in the selected etiology,
            or <b>All Signatures</b> to see signatures for every etiology.
          </p>
        </Col>
      </Row>
      <Row className="justify-content-center mb-3">
        <Col sm="auto">
          <Check id="selectedEtiology">
            <Check.Input
              type="radio"
              checked={!showAllSignatures}
              onChange={() => mergeEtiology({ showAllSignatures: false })}
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
              checked={showAllSignatures}
              onChange={() => mergeEtiology({ showAllSignatures: true })}
            />
            <Check.Label className="font-weight-normal">
              All Signatures
            </Check.Label>
          </Check>
        </Col>
      </Row>

      {isFetching && <LoadingOverlay active="true" />}
      {error && (
        <p>
          There was an error retrieving signature options. Please try again at a
          later time.
        </p>
      )}
      {thumbnails && (
        <>
          <Row
            className={`justify-content-center ${
              showAllSignatures ? 'd-none' : ''
            }`}
          >
            {unqiueSignatures
              .filter((e) => e.etiology == etiology)
              .sort(naturalSort)
              .sort(profileSort)
              .map((e, i) => {
                return (
                  <Col key={i} md="2" sm="4" className="mb-3">
                    <div
                      className={`sigIcon border rounded ${
                        e.signature == signature ? 'active' : ''
                      }`}
                      title={`${etiology} - ${e.signature}`}
                      onClick={() =>
                        mergeEtiology({
                          signature: e.signature,
                        })
                      }
                    >
                      <img
                        src={thumbnails[e.signature.replaceAll(/\//g, ':')]}
                        className="w-100"
                        // height="110"
                        alt=""
                      />
                      <strong className="sigLabel">{e.signature}</strong>
                    </div>
                  </Col>
                );
              })}
          </Row>
          <Row
            className={`px-2 justify-content-center ${
              !showAllSignatures ? 'd-none' : ''
            }`}
          >
            {unqiueSignatures
              .sort(naturalSort)
              .sort(profileSort)
              .map((e, i) => (
                <Col key={i} lg="1" md="3" sm="4" className="mb-2 px-1">
                  <div
                    onClick={() =>
                      mergeEtiology({
                        etiology: e.etiology,
                        signature: e.signature,
                      })
                    }
                    className={`sigIcon border rounded ${
                      etiology != e.etiology
                        ? 'inactive'
                        : signature == e.signature
                        ? 'active'
                        : ''
                    }`}
                    title={`${e.etiology} - ${e.signature}`}
                  >
                    <img
                      src={thumbnails[e.signature.replaceAll(/\//g, ':')]}
                      className="w-100"
                      // height="70"
                      alt=""
                    />
                    <strong className="sigLabel" style={{ fontSize: '0.8rem' }}>
                      {e.signature}
                    </strong>
                  </div>
                </Col>
              ))}
          </Row>
        </>
      )}
      {refSigThumbnails && category == 'CancerSpecificSignatures' && (
        <div className="my-3">
          <h5 className="separator">Reference Signatures</h5>
          {!signature && (
            <p className="text-center text-muted">Select a Signature</p>
          )}
          <Row className="justify-content-center">
            {data
              .filter((e) => e.signature == signature)
              .map((e, i) => {
                return (
                  <Col key={i} md="2" sm="4" className="mb-3">
                    <div
                      className={`sigIcon border rounded ${
                        referenceSignature == e.json.referenceSignature
                          ? 'active'
                          : ''
                      }`}
                      title={`${e.signature} - ${e.json.referenceSignature}`}
                      onClick={() =>
                        mergeEtiology({
                          referenceSignature: e.json.referenceSignature,
                        })
                      }
                    >
                      <img
                        src={
                          refSigThumbnails[
                            e.json.referenceSignature.replaceAll(/\//g, ':')
                          ]
                        }
                        className="w-100"
                        // height="110"
                        alt=""
                      />
                      <strong className="sigLabel">
                        {`${e.json.referenceSignature} (${e.json.refSigProportion})`}
                      </strong>
                    </div>
                  </Col>
                );
              })}
          </Row>
        </div>
      )}
    </div>
  );
}
