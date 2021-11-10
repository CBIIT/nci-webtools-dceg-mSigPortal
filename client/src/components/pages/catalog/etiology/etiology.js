import React, { useEffect } from 'react';
import { Row, Col, Button, Form, Card } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';
import { getJSON } from '../../../../services/utils';
import Plot from '../../../controls/plot/plot';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import './etiology.scss';

const actions = { ...catalogActions, ...modalActions };
const { Check } = Form;

export default function Etiology() {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeCatalog({ etiology: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    category,
    etiology,
    selectedSignature,
    tissue,
    refSig,
    study,
    all,
    data,
    thumbnails,
    tissueThumbnails,
    refSigThumbnails,
    profileURL,
    exposureURL,
    strandbiasURL,
    tissueURL,
    refSigURL,
    loading,
  } = useSelector((state) => state.catalog.etiology);

  const categories = [
    {
      name: 'Cosmic Mutational Signatures (v3.2)',
      file: 'Etiology_cosmic.json',
    },
    {
      name: 'Environmental Mutagenesis',
      file: 'Etiology_enviromental_mutagenesis.json',
    },
    { name: 'Gene Edits', file: 'Etiology_gene_edits.json' },
    {
      name: 'Cancer Specific Signature',
      file: 'Etiology_cancer_specific_signatures.json',
    },
    { name: 'Others', file: 'Etiology_others.json' },
  ];

  useEffect(() => {
    dispatch(actions.mergeCatalog({ catalog: { displayTab: 'etiology' } }));
  }, []);

  useEffect(() => {
    const getData = async () => {
      mergeEtiology({ loading: true });
      try {
        const file = categories.filter(({ name }) => name == category)[0].file;
        const data = await getJSON(`Etiology/${file}`);

        const etiologyOptions = [
          ...new Set(data.map(({ Etiology }) => Etiology)),
        ];
        const studyOptions = [...new Set(data.map(({ Study }) => Study))];

        mergeEtiology({
          data: { [category]: data },
          etiology: etiologyOptions[0],
          study: studyOptions[0],
        });
      } catch (e) {
        mergeError(e.message);
        console.error(e);
      }
      mergeEtiology({ loading: false });
    };
    if (!data[category]) {
      if (!loading) getData();
    } else {
      const etiologyOptions = [
        ...new Set(data[category].map(({ Etiology }) => Etiology)),
      ];
      const studyOptions = [
        ...new Set(data[category].map(({ Study }) => Study)),
      ];

      mergeEtiology({
        etiology: etiologyOptions[0],
        study: studyOptions[0],
      });
    }
  }, [category, loading]);

  // check if plot exists
  useEffect(() => {
    const getImageS3 = (path) =>
      fetch(`api/getImageS3`, {
        method: 'POST',
        headers: {
          Accept: 'image/svg',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          path: `msigportal/Database/Etiology/${path}`,
        }),
      }).then(async (res) => {
        if (res.ok) {
          const blob = await res.blob();
          return URL.createObjectURL(blob);
        } else {
          return '';
        }
      });

    const getPlots = async () => {
      if (profileURL) URL.revokeObjectURL(profileURL);
      if (exposureURL) URL.revokeObjectURL(exposureURL);
      if (strandbiasURL) URL.revokeObjectURL(strandbiasURL);

      const [sig, tmb, strandBias] = await Promise.all([
        getImageS3(`Profile/${fixFile(selectedSignature)}.svg`),
        getImageS3(`Exposure/${fixFile(`${selectedSignature}_${study}`)}.svg`),
        getImageS3(`Profile_StrandBias/${fixFile(selectedSignature)}.svg`),
      ]);
      mergeEtiology({
        profileURL: sig,
        exposureURL: tmb,
        strandbiasURL: strandBias,
      });
    };

    const getCancerSpecificPlots = async () => {
      if (tissueURL) URL.revokeObjectURL(tissueURL);
      if (refSigURL) URL.revokeObjectURL(refSigURL);
      if (exposureURL) URL.revokeObjectURL(exposureURL);
      const [tissuePlot, refSigPlot, exposurePlot] = await Promise.all([
        getImageS3(`Profile/${fixFile(tissue)}.svg`),
        getImageS3(`Profile/${fixFile(refSig)}.svg`),
        getImageS3(`Exposure/${fixFile(tissue)}_PCAWG.svg`),
      ]);
      mergeEtiology({
        tissueURL: tissuePlot,
        refSigURL: refSigPlot,
        exposureURL: exposurePlot,
      });
    };

    if (selectedSignature) getPlots();
    else if (tissue && refSig) getCancerSpecificPlots();
  }, [selectedSignature, study, tissue, refSig]);

  // get thumbnails for standard categories
  useEffect(() => {
    if (
      data[category] &&
      data[category].length &&
      !thumbnails[category] &&
      category != 'Cancer Specific Signature' &&
      !loading
    ) {
      const signatures = data[category]
        .slice()
        .filter(
          (value, index, array) =>
            array.findIndex(
              (t) => t['Signature Name'] === value['Signature Name']
            ) === index ||
            array.findIndex((t) => t['Signature'] === value['Signature']) ===
              index
        )
        .sort(naturalSort)
        .sort(profileSort);

      getThumbnails(signatures);
    }
  }, [data, category, loading]);

  // get tissue thumbnails
  useEffect(() => {
    if (
      data[category] &&
      data[category].length &&
      !tissueThumbnails.length &&
      category == 'Cancer Specific Signature'
    ) {
      const uniqueTissues = Object.values(
        data[category].reduce((c, e) => {
          if (!c[e['Tissue Specific Signature']])
            c[e['Tissue Specific Signature']] = e;
          return c;
        }, {})
      ).sort(naturalSort);

      getTissueThumbnails(uniqueTissues);
    }
  }, [data, category]);

  // get refsig thumbnails
  useEffect(() => {
    if (
      data[category] &&
      data[category].length &&
      !refSigThumbnails.length &&
      category == 'Cancer Specific Signature'
    ) {
      const refsig = data[category].slice().sort(naturalSort);

      getRefSigThumbnails(refsig);
    }
  }, [data]);

  function profileSort(a, b) {
    const sigOrder = [/SBS/g, /DBS/g, /ID/g];
    let c = 0,
      d = 0;

    sigOrder.forEach((profile, i) => {
      const sigA = a['Signature Name'] || a.Signature;
      const sigB = b['Signature Name'] || b.Signature;
      if (sigA.match(profile)) c = i;
      if (sigB.match(profile)) d = i;
    });

    return c - d;
  }

  // replace forward slash with colon for proper fs traversal
  function fixFile(filename) {
    return filename.replace(/\//g, ':');
  }

  function naturalSort(a, b) {
    const sigA = a['Signature Name'] || a.Signature || false;
    const sigB = b['Signature Name'] || b.Signature || false;
    if (sigA)
      return sigA.localeCompare(sigB, undefined, {
        numeric: true,
        sensitivity: 'base',
      });
    else if (a['Tissue Specific Signature'])
      return a['Tissue Specific Signature'].localeCompare(
        b['Tissue Specific Signature'],
        undefined,
        {
          numeric: true,
          sensitivity: 'base',
        }
      );
    else {
      return 0;
    }
  }

  function getCategories() {
    return (
      <Row className="justify-content-center mb-3">
        {categories.map(({ name, file }) => (
          <Col key={name} lg="2" md="3" sm="4" className="mb-3 d-flex">
            <Button
              size="sm"
              variant="dark"
              onClick={
                file && name != category
                  ? async () => {
                      mergeEtiology({
                        category: name,
                        etiology: '',
                        selectedSignature: '',
                        study: '',
                      });
                    }
                  : () => {}
              }
              className={category != name ? 'disabled' : ''}
              block
            >
              {name}
            </Button>
          </Col>
        ))}
      </Row>
    );
  }

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
                  lg={category == 'Gene Edits' ? '1' : '2'}
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
                        selectedSignature: '',
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

  async function getImageBatch(keyArr) {
    return (
      await fetch(`api/getImageS3Batch`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          keys: keyArr,
        }),
      })
    ).json();
  }

  async function getThumbnails(signatures) {
    try {
      const keys = signatures.map(
        (signature) =>
          `msigportal/Database/Etiology/Profile_logo/${fixFile(
            signature['Signature Name'] || signature.Signature
          )}.svg`
      );

      const images = await getImageBatch(keys);

      const thumbnails = signatures.map((signature, i) => ({
        Etiology: signature.Etiology,
        Study: signature.Study,
        signatureName: signature['Signature Name'] || signature.Signature,
        thumbnailURL: images[i],
      }));

      mergeEtiology({ thumbnails: { [category]: thumbnails } });
    } catch (err) {
      mergeError(err.message);
      console.error(err);
    }
  }

  async function getTissueThumbnails(tissues) {
    try {
      const keys = tissues.map(
        (tissue) =>
          `msigportal/Database/Etiology/Profile_logo/${fixFile(
            tissue['Tissue Specific Signature']
          )}.svg`
      );

      const images = await getImageBatch(keys);

      const thumbnails = tissues.map((tissue, i) => ({
        Etiology: tissue.Etiology,
        'Tissue Specific Signature': tissue['Tissue Specific Signature'],
        'Ref Signature': tissue['Ref Signature'],
        thumbnailURL: images[i],
      }));

      mergeEtiology({ tissueThumbnails: thumbnails });
    } catch (err) {
      mergeError(err.message);
      console.error(err);
    }
  }

  async function getRefSigThumbnails(refSigs) {
    try {
      const keys = refSigs.map(
        (refSig) =>
          `msigportal/Database/Etiology/Profile_logo/${fixFile(
            refSig['Ref Signature']
          )}.svg`
      );

      const images = await getImageBatch(keys);

      const thumbnails = refSigs.map((refSig, i) => ({
        Etiology: refSig.Etiology,
        'Tissue Specific Signature': refSig['Tissue Specific Signature'],
        'Ref Signature': refSig['Ref Signature'],
        'RefSig Proportion': refSig['RefSig Proportion'],
        thumbnailURL: images[i],
      }));

      mergeEtiology({ refSigThumbnails: thumbnails });
    } catch (err) {
      mergeError(err.message);
      console.error(err);
    }
  }

  function standardView() {
    function getAllSignatures() {
      if (thumbnails[category] && thumbnails[category].length) {
        return thumbnails[category].map(
          ({ Etiology, signatureName, thumbnailURL, Source_URL }, index) => (
            <Col key={index} lg="1" md="3" sm="4" className="mb-2 px-1">
              <div
                onClick={() =>
                  mergeEtiology({
                    etiology: Etiology,
                    selectedSignature: signatureName,
                  })
                }
                className={`sigIcon border rounded ${
                  etiology != Etiology
                    ? 'inactive'
                    : signatureName == selectedSignature
                    ? 'active'
                    : ''
                }`}
                title={`${Etiology} - ${signatureName}`}
              >
                <img
                  src={thumbnailURL}
                  className="w-100"
                  // height="70"
                  alt={signatureName}
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
                    signatureName == selectedSignature ? 'active' : ''
                  }`}
                  title={`${etiology} - ${signatureName}`}
                  onClick={() =>
                    mergeEtiology({
                      selectedSignature: signatureName,
                    })
                  }
                >
                  <img
                    src={thumbnailURL}
                    className="w-100"
                    // height="110"
                    alt={signatureName}
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
                  Etiology == etiology && signatureName == selectedSignature
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
      if (data[category] && data[category].length && selectedSignature) {
        let info = data[category].filter(
          (signature) =>
            signature['Signature Name'] == selectedSignature ||
            signature['Signature'] == selectedSignature
        );
        if (info.length) {
          info = info[0];

          return (
            <div>
              <div>
                <strong>Signature Name: </strong>
                {info['Signature Name']}
              </div>
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
              {info['Cell Line'] && (
                <div>
                  <strong>Cell Line: </strong>
                  {info['Cell Line']}
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
                  <Plot
                    className="p-3 border rounded mb-3"
                    height={'500px'}
                    plotPath={profileURL}
                    cacheBreaker={false}
                  />

                  {info.Description_strandbias && (
                    <p>{info.Description_strandbias}</p>
                  )}
                  {info.Description_strandbias && strandbiasURL && (
                    <Plot
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
                          megabase (log10) attributed to each mutational
                          signature in samples where the signature is present.
                          Only those cancer types with tumors in which signature
                          activity is attributed are shown. The numbers below
                          the dots for each cancer type indicate the number of
                          tumors in which the signatures was attributed (above
                          the horizontal bar, in blue) and the total number of
                          tumors analyzed (below the blue bar, in green).
                        </p>
                        <Row className="justify-content-center">
                          {getStudy()}
                        </Row>
                        {exposureURL.length > 0 && (
                          <Plot
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
              Error: No data found for {selectedSignature}
            </p>
          );
        }
      }
    }

    return (
      <div>
        <div className="mb-3">
          <h5 className="separator">Etiologies</h5>
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
          <Row
            className={`px-2 justify-content-center ${!all ? 'd-none' : ''}`}
          >
            {getAllSignatures()}
          </Row>
          <div className="p-3">{getInfo()}</div>
        </div>
      </div>
    );
  }

  function cancerSpecificView() {
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
                  alt={tissueSig['Tissue Specific Signature']}
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
                    alt={tissueSig['Tissue Specific Signature']}
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
                        alt={v['Ref Signature']}
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
          <Row
            className={`px-2 justify-content-center ${!all ? 'd-none' : ''}`}
          >
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

  return (
    <div className="p-4 bg-white border rounded">
      <h5 className="separator">Categories</h5>
      {getCategories()}
      {category != 'Cancer Specific Signature' && standardView()}
      {category == 'Cancer Specific Signature' && cancerSpecificView()}
    </div>
  );
}
