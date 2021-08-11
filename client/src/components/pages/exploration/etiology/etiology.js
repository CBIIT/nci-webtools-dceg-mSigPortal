import React, { useEffect } from 'react';
import { Row, Col, Button, Form } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as etiologyActions } from '../../../../services/store/etiology';
import { actions as explorationActions } from '../../../../services/store/exploration';
import { actions as modalActions } from '../../../../services/store/modal';
import { getJSON } from '../../../../services/utils';
import Plot from '../../../controls/plot/plot';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import './etiology.scss';

const actions = { ...etiologyActions, ...explorationActions, ...modalActions };
const { Group, Check } = Form;

export default function Etiology() {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeEtiology({ etiologyState: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    category,
    aetiology,
    signature,
    tissue,
    refSig,
    study,
    all,
    data,
    thumbnails,
    tissueThumbnails,
    refSigThumbnails,
    selectedSource,
    profileURL,
    exposureURL,
    tissueURL,
    refSigURL,
  } = useSelector((state) => state.etiology.etiologyState);

  const categories = [
    { name: 'Cosmic Mutational Signatures', file: 'Aetiology_cosmic.json' },
    {
      name: 'Environmental Mutagenesis',
      file: 'Aetiology_enviromental_mutagenesis.json',
    },
    { name: 'Gene Edits', file: 'Aetiology_gene_edits.json' },
    {
      name: 'Cancer Specific Signature',
      file: 'Aetiology_cancer_specific_signatures.json',
    },
    { name: 'Others', file: '' },
  ];

  useEffect(() => {
    dispatch(
      actions.mergeExploration({ exploration: { displayTab: 'etiology' } })
    );
  }, []);

  useEffect(() => {
    const getData = async () => {
      try {
        const file = categories.filter(({ name }) => name == category)[0].file;
        const data = await getJSON(`Aetiology/${file}`);

        const aetiologyOptions = [
          ...new Set(data.map(({ Aetiology }) => Aetiology)),
        ];
        const studyOptions = [...new Set(data.map(({ Study }) => Study))];

        mergeEtiology({
          data: { [category]: data },
          aetiology: aetiologyOptions[0],
          study: studyOptions[0],
        });
      } catch (e) {
        mergeError(e);
      }
    };
    if (!data[category]) {
      getData();
    } else {
      const aetiologyOptions = [
        ...new Set(data[category].map(({ Aetiology }) => Aetiology)),
      ];
      const studyOptions = [
        ...new Set(data[category].map(({ Study }) => Study)),
      ];

      mergeEtiology({
        aetiology: aetiologyOptions[0],
        study: studyOptions[0],
      });
    }
  }, [category]);

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
          path: `msigportal/Database/Aetiology/${path}`,
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
      const [sig, tmb] = await Promise.all([
        getImageS3(`Profile/${fixFile(signature)}.svg`),
        getImageS3(`Exposure/${fixFile(`${signature}_${study}`)}.svg`),
      ]);
      mergeEtiology({ profileURL: sig, exposureURL: tmb });
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

    if (signature && study) getPlots();
    else if (tissue && refSig) getCancerSpecificPlots();
  }, [signature, study, tissue, refSig]);

  // get thumbnails for standard categories
  useEffect(() => {
    if (
      data[category] &&
      data[category].length &&
      !thumbnails[category] &&
      category != 'Cancer Specific Signature'
    ) {
      const signatures = data[category]
        .slice()
        .filter(({ Study }) => Study == study)
        .sort(naturalSort)
        .sort(profileSort);

      getThumbnails(signatures);
    }
  }, [data]);

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
  }, [data]);

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
    const sigOrder = ['SBS', 'DBS', 'ID'];
    let c = 0,
      d = 0;

    sigOrder.forEach((profile, i) => {
      if (a.Signature && a.Signature.includes(profile)) {
        c = i;
      }
      if (b.Signature && b.Signature.includes(profile)) {
        d = i;
      }
    });

    return c - d;
  }

  // replace forward slash with colon for proper fs traversal
  function fixFile(filename) {
    return filename.replace(/\//g, ':');
  }

  function naturalSort(a, b) {
    if (a.Signature)
      return a.Signature.localeCompare(b.Signature, undefined, {
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
    return categories.map(({ name, file }) => (
      <Col key={name} lg="2" md="3" sm="4" className="mb-3 d-flex">
        <Button
          size="sm"
          variant="dark"
          onClick={
            file && name != category
              ? async () => {
                  mergeEtiology({
                    category: name,
                    aetiology: '',
                    signature: '',
                    study: '',
                    selectedSource: '',
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
    ));
  }

  function getAetiologies() {
    if (data[category] && data[category].length) {
      return (
        <Row className="justify-content-center">
          {[...new Set(data[category].map((obj) => obj.Aetiology))]
            .sort()
            .map((Aetiology) => (
              <Col
                key={Aetiology}
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
                      aetiology: Aetiology,
                      signature: '',
                      selectedSource: '',
                    })
                  }
                  className={aetiology != Aetiology ? 'disabled' : ''}
                  block
                >
                  {Aetiology}
                </Button>
              </Col>
            ))}
        </Row>
      );
    } else {
      return (
        <div className="my-5">
          <LoadingOverlay active={true} />
        </div>
      );
    }
  }
  function getCancerAetiology() {
    if (data[category] && data[category].length) {
      return (
        <Row className="justify-content-center">
          {[...new Set(data[category].map((obj) => obj.Aetiology))]
            .sort()
            .map((Aetiology) => (
              <Col key={Aetiology} lg="2" md="3" sm="4" className="mb-3 d-flex">
                <Button
                  size="sm"
                  variant="dark"
                  onClick={() =>
                    mergeEtiology({
                      aetiology: Aetiology,
                      tissue: '',
                      refSig: '',
                      selectedSource: '',
                    })
                  }
                  className={aetiology != Aetiology ? 'disabled' : ''}
                  block
                >
                  {Aetiology}
                </Button>
              </Col>
            ))}
        </Row>
      );
    } else {
      return (
        <div className="my-5">
          <LoadingOverlay active={true} />
        </div>
      );
    }
  }

  async function getThumbnails(signatures) {
    try {
      const newThumbnails = await Promise.all(
        signatures.map(
          ({ Aetiology, Study, Signature, Source_URL, URL: alias_URL }) =>
            fetch(`api/getImageS3`, {
              method: 'POST',
              headers: {
                Accept: 'image/svg',
                'Content-Type': 'application/json',
              },
              body: JSON.stringify({
                path: `msigportal/Database/Aetiology/Profile_logo/${fixFile(
                  Signature
                )}.svg`,
              }),
            }).then(async (res) => {
              const blob = await res.blob();
              return {
                Aetiology: Aetiology,
                Study: Study,
                Signature: Signature,
                // use the Source_URL as a unique identifier
                Source_URL: Source_URL || alias_URL,
                url: URL.createObjectURL(blob),
              };
            })
        )
      );

      mergeEtiology({ thumbnails: { [category]: newThumbnails } });
    } catch (err) {
      mergeError(err.message);
    }
  }

  async function getTissueThumbnails(tissues) {
    try {
      tissueThumbnails.forEach((t) => URL.revokeObjectURL(t));

      const newThumbnails = await Promise.all(
        tissues.map((tissue) =>
          fetch(`api/getImageS3`, {
            method: 'POST',
            headers: {
              Accept: 'image/svg',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              path: `msigportal/Database/Aetiology/Profile_logo/${fixFile(
                tissue['Tissue Specific Signature']
              )}.svg`,
            }),
          }).then(async (res) => {
            const blob = await res.blob();
            return {
              Aetiology: tissue.Aetiology,
              'Tissue Specific Signature': tissue['Tissue Specific Signature'],
              'Ref Signature': tissue['Ref Signature'],
              url: URL.createObjectURL(blob),
            };
          })
        )
      );
      mergeEtiology({ tissueThumbnails: newThumbnails });
    } catch (err) {
      mergeError(err.message);
    }
  }

  async function getRefSigThumbnails(refSigs) {
    try {
      refSigThumbnails.forEach((t) => URL.revokeObjectURL(t));

      const newThumbnails = await Promise.all(
        refSigs.map((refSig) =>
          fetch(`api/getImageS3`, {
            method: 'POST',
            headers: {
              Accept: 'image/svg',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              path: `msigportal/Database/Aetiology/Profile_logo/${fixFile(
                refSig['Ref Signature']
              )}.svg`,
            }),
          }).then(async (res) => {
            const blob = await res.blob();
            return {
              Aetiology: refSig.Aetiology,
              'Tissue Specific Signature': refSig['Tissue Specific Signature'],
              'Ref Signature': refSig['Ref Signature'],
              'RefSig Proportion': refSig['RefSig Proportion'],
              url: URL.createObjectURL(blob),
            };
          })
        )
      );
      mergeEtiology({ refSigThumbnails: newThumbnails });
    } catch (err) {
      mergeError(err.message);
    }
  }

  function standardView() {
    function getAllSignatures() {
      if (thumbnails[category] && thumbnails[category].length) {
        return thumbnails[category].map(
          ({ Aetiology, Signature, url, Source_URL }, index) => (
            <Col key={index} lg="1" md="3" sm="4" className="mb-2 px-1">
              <div
                onClick={() =>
                  mergeEtiology({
                    aetiology: Aetiology,
                    signature: Signature,
                    selectedSource: Source_URL,
                  })
                }
                className={`sigIcon border rounded ${
                  aetiology != Aetiology
                    ? 'inactive'
                    : selectedSource == Source_URL
                    ? 'active'
                    : ''
                }`}
                title={`${Aetiology} - ${Signature}`}
              >
                <img
                  src={url}
                  className="w-100"
                  // height="70"
                  alt={Signature}
                />
                <div className="sigLabel">
                  <strong style={{ fontSize: '0.8rem' }}>{Signature}</strong>
                </div>
              </div>
            </Col>
          )
        );
      } else {
        return (
          <div className="my-5">
            <LoadingOverlay active={true} />
          </div>
        );
      }
    }

    function getSignatures() {
      if (thumbnails[category] && thumbnails[category].length) {
        return thumbnails[category]
          .filter(({ Aetiology }) => Aetiology == aetiology)
          .map(({ Signature, url, Source_URL }, index) => {
            return (
              <Col key={index} md="2" sm="4" className="mb-3">
                <div
                  className={`sigIcon border rounded ${
                    selectedSource == Source_URL ? 'active' : ''
                  }`}
                  title={`${aetiology} - ${Signature}`}
                  onClick={() =>
                    mergeEtiology({
                      signature: Signature,
                      selectedSource: Source_URL,
                    })
                  }
                >
                  <img
                    src={url}
                    className="w-100"
                    // height="110"
                    alt={Signature}
                  />
                  <div className="sigLabel">
                    <strong className="sigLabel">{Signature}</strong>
                  </div>
                </div>
              </Col>
            );
          });
      } else {
        return (
          <div className="my-5">
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
                ({ Aetiology, Signature }) =>
                  Aetiology == aetiology && Signature == signature
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
      if (data[category] && data[category].length && selectedSource) {
        let info = data[category].filter(
          ({ Source_URL, URL: alias_URL }) =>
            Source_URL == selectedSource || alias_URL == selectedSource
        );
        if (info.length) {
          info = info[0];

          return (
            <div>
              <div>
                <strong>Signature Name: </strong>
                {signature}
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
              {info.Source && (
                <div>
                  <strong>Source: </strong>
                  <a
                    href={info.URL || info.Source_URL}
                    target="_blank"
                    rel="noreferrer"
                  >
                    {info.Source}
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
                  />

                  {category == 'Cosmic Mutational Signatures' && (
                    <>
                      <>
                        <Row className="justify-content-center">
                          {getStudy()}
                        </Row>
                        <p>
                          Select the cancer study to review the TMB of selected
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
                        {exposureURL.length ? (
                          <Plot
                            className="p-3 border"
                            height={'600px'}
                            plotPath={exposureURL}
                          />
                        ) : (
                          <div className="p-3 border">
                            <p>
                              A signature was not detected in any sample of the
                              selected study
                            </p>
                          </div>
                        )}
                      </>
                    </>
                  )}
                </div>
              ) : (
                <div className="my-5">
                  <LoadingOverlay active={true} />
                </div>
              )}
            </div>
          );
        } else {
          return (
            <p className="d-flex justify-content-center text-muted">
              Error: No data found for {selectedSource}
            </p>
          );
        }
      } else {
        return (
          <p className="d-flex justify-content-center text-muted">
            Select a category, etiology, and signature to view more info
          </p>
        );
      }
    }

    return (
      <div>
        <div className="mx-auto p-3 pt-0">
          <Row>
            <Col sm="4">
              <strong>Etiologies</strong>
            </Col>
            <Col sm="4">
              <Group className="d-flex justify-content-center">
                <Check inline id="selectedAetiology">
                  <Check.Input
                    type="radio"
                    checked={all == false}
                    onChange={() => mergeEtiology({ all: false })}
                  />
                  <Check.Label className="font-weight-normal">
                    Selected Etiology
                  </Check.Label>
                </Check>
                <Check inline id="allAetiologies">
                  <Check.Input
                    type="radio"
                    checked={all == true}
                    onChange={() => mergeEtiology({ all: true })}
                  />
                  <Check.Label className="font-weight-normal">
                    All Etiologies
                  </Check.Label>
                </Check>
              </Group>
            </Col>
          </Row>
          {getAetiologies()}
        </div>
        <hr />
        <div className="mx-auto p-3 pt-0">
          <strong>Signatures</strong>
          <Row className={`justify-content-center ${all ? 'd-none' : ''}`}>
            {getSignatures()}
          </Row>
          <Row
            className={`px-2 justify-content-center ${!all ? 'd-none' : ''}`}
          >
            {getAllSignatures()}
          </Row>
        </div>
        <hr />
        <div className="mx-auto p-3">{getInfo()}</div>
      </div>
    );
  }

  function cancerSpecificView() {
    function getAllTissues() {
      if (tissueThumbnails.length) {
        return tissueThumbnails.map((tissue, index) => {
          return (
            <Col key={index} lg="1" md="2" sm="4" className="mb-2 px-1">
              <div
                className={`sigIcon border rounded ${
                  aetiology != tissue.Aetiology
                    ? 'inactive'
                    : selectedSource ==
                      tissue['Tissue Specific Signature'] +
                        tissue['Ref Signature']
                    ? 'active'
                    : ''
                }`}
                title={`${aetiology} - ${tissue['Tissue Specific Signature']}`}
                onClick={() =>
                  mergeEtiology({
                    aetiology: tissue.Aetiology,
                    tissue: tissue['Tissue Specific Signature'],
                    selectedSource:
                      tissue['Tissue Specific Signature'] +
                      tissue['Ref Signature'],
                  })
                }
              >
                <img
                  src={tissue.url}
                  className="w-100"
                  // height="110"
                  alt={tissue['Tissue Specific Signature']}
                />
                <div className="sigLabel">
                  <strong className="sigLabel">
                    {tissue['Tissue Specific Signature']}
                  </strong>
                </div>
              </div>
            </Col>
          );
        });
      } else {
        return (
          <div className="my-5">
            <LoadingOverlay active={true} />
          </div>
        );
      }
    }

    function getTissues() {
      if (tissueThumbnails.length) {
        return tissueThumbnails
          .filter(({ Aetiology }) => Aetiology == aetiology)
          .map((tissue, index) => {
            return (
              <Col key={index} md="2" sm="4" className="mb-3">
                <div
                  className={`sigIcon border rounded ${
                    selectedSource ==
                    tissue['Tissue Specific Signature'] +
                      tissue['Ref Signature']
                      ? 'active'
                      : ''
                  }`}
                  title={`${aetiology} - ${tissue['Tissue Specific Signature']}`}
                  onClick={() =>
                    mergeEtiology({
                      tissue: tissue['Tissue Specific Signature'],
                      refSig: '',
                      selectedSource:
                        tissue['Tissue Specific Signature'] +
                        tissue['Ref Signature'],
                    })
                  }
                >
                  <img
                    src={tissue.url}
                    className="w-100"
                    // height="110"
                    alt={tissue['Tissue Specific Signature']}
                  />
                  <div className="sigLabel">
                    <strong className="sigLabel">
                      {tissue['Tissue Specific Signature']}
                    </strong>
                  </div>
                </div>
              </Col>
            );
          });
      } else {
        return (
          <div className="my-5">
            <LoadingOverlay active={true} />
          </div>
        );
      }
    }

    function getRefSig() {
      if (refSigThumbnails.length && tissue) {
        return refSigThumbnails
          .filter(
            (v) =>
              v.Aetiology == aetiology &&
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
                  onClick={() => mergeEtiology({ refSig: v['Ref Signature'] })}
                >
                  <img
                    src={v.url}
                    className="w-100"
                    // height="110"
                    alt={v['Ref Signature']}
                  />
                  <div className="sigLabel">
                    <strong className="sigLabel">
                      {`${v['Ref Signature']} (${v['RefSig Proportion']})`}
                    </strong>
                  </div>
                </div>
              </Col>
            );
          });
      } else {
        return <p className="text-muted">Select a signature</p>;
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
      } else {
        return (
          <p className="d-flex justify-content-center text-muted">
            Select a Reference Signature
          </p>
        );
      }
    }
    return (
      <div>
        <div className="mx-auto p-3">
          <Row>
            <Col sm="4">
              <strong>Etiologies</strong>
            </Col>
            <Col sm="4">
              <Group className="d-flex justify-content-center">
                <Check inline id="selectedAetiology">
                  <Check.Input
                    type="radio"
                    checked={all == false}
                    onChange={() => mergeEtiology({ all: false })}
                  />
                  <Check.Label className="font-weight-normal">
                    Selected Etiology
                  </Check.Label>
                </Check>
                <Check inline id="allAetiologies">
                  <Check.Input
                    type="radio"
                    checked={all == true}
                    onChange={() => mergeEtiology({ all: true })}
                  />
                  <Check.Label className="font-weight-normal">
                    All Etiologies
                  </Check.Label>
                </Check>
              </Group>
            </Col>
          </Row>
          {getCancerAetiology()}
        </div>
        <hr />
        <div className="mx-auto p-3">
          <strong>Tissue Specific Signatures</strong>
          <Row className={`justify-content-center ${all ? 'd-none' : ''}`}>
            {getTissues()}
          </Row>
          <Row
            className={`px-2 justify-content-center ${!all ? 'd-none' : ''}`}
          >
            {getAllTissues()}
          </Row>
        </div>
        <hr />
        <div className="mx-auto p-3">
          <strong>Reference Signatures</strong>
          <Row className="justify-content-center">{getRefSig()}</Row>
        </div>
        <hr />
        <div className="mx-auto p-3">{getCancerSpecificInfo()}</div>
      </div>
    );
  }

  return (
    <div className="bg-white border rounded">
      <div className="mx-auto p-3">
        <strong className="mb-3">Categories</strong>
        <Row className="justify-content-center">{getCategories()}</Row>
      </div>
      <hr />
      {category != 'Cancer Specific Signature' && standardView()}
      {category == 'Cancer Specific Signature' && cancerSpecificView()}
    </div>
  );
}
