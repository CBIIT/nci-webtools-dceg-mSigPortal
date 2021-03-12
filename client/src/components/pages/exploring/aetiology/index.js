import React, { useEffect } from 'react';
import { Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import './aetiology.scss';

const actions = { ...exploringActions, ...modalActions };

export default function Aetiology() {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const mergeAetiology = (state) =>
    dispatch(actions.mergeExploring({ aetiology: state }));
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
  } = exploring.aetiology;

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
    const getData = async () => {
      try {
        const file = categories.filter(({ name }) => name == category)[0].file;
        const data = await (
          await fetch(`api/getFileS3`, {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              path: `msigportal/Database/Aetiology/${file}`,
            }),
          })
        ).json();

        const aetiologyOptions = [
          ...new Set(data.map(({ Aetiology }) => Aetiology)),
        ];
        const studyOptions = [...new Set(data.map(({ Study }) => Study))];

        mergeAetiology({
          data: data,
          aetiology: aetiologyOptions[0],
          study: studyOptions[0],
        });
      } catch (_) {
        mergeError('Could not fetch Aetiology Info');
      }
    };
    if (!data.length) getData();
  }, [data]);

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
          path: path,
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
        getImageS3(`msigportal/Database/Aetiology/Profile/${signature}.svg`),
        getImageS3(
          `msigportal/Database/Aetiology/Exposure/${signature}_${study}.svg`
        ),
      ]);
      mergeAetiology({ profileURL: sig, exposureURL: tmb });
    };

    const getCancerSpecificPlots = async () => {
      if (tissueURL) URL.revokeObjectURL(tissueURL);
      if (refSigURL) URL.revokeObjectURL(refSigURL);
      if (exposureURL) URL.revokeObjectURL(exposureURL);
      const [tissuePlot, refSigPlot, exposurePlot] = await Promise.all([
        getImageS3(`msigportal/Database/Aetiology/Profile/${tissue}.svg`),
        getImageS3(`msigportal/Database/Aetiology/Profile/${refSig}.svg`),
        getImageS3(
          `msigportal/Database/Aetiology/Exposure/${tissue}_PCAWG.svg`
        ),
      ]);
      mergeAetiology({
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
    if (data.length && category != 'Cancer Specific Signature') {
      const signatures = data
        .slice()
        .filter(({ Study }) => Study == study)
        .sort(naturalSort)
        .sort(profileSort);

      getThumbnails(signatures);
    }
  }, [data]);

  // get tissue thumbnails
  useEffect(() => {
    if (data.length && category == 'Cancer Specific Signature') {
      const uniqueTissues = Object.values(
        data.reduce((c, e) => {
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
    if (data.length && category == 'Cancer Specific Signature') {
      const refsig = data.slice().sort(naturalSort);

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
          variant="primary"
          onClick={
            file && name != category
              ? async () => {
                  mergeAetiology({
                    category: name,
                    aetiology: '',
                    signature: '',
                    study: '',
                    selectedSource: '',
                    data: [],
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
    if (data.length) {
      const aetiologies = [...new Set(data.map((obj) => obj.Aetiology))]
        .sort()
        .map((Aetiology) => (
          <Col key={Aetiology} lg="2" md="3" sm="4" className="mb-3 d-flex">
            <Button
              size="sm"
              variant="dark"
              onClick={() =>
                mergeAetiology({
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
        ));

      return (
        <>
          <Row className="justify-content-center">{aetiologies}</Row>
          <Row className="justify-content-center">
            <Col lg="2" md="3" sm="4" className="mb-3 d-flex">
              <Button
                size="sm"
                variant="primary"
                onClick={() => mergeAetiology({ all: false })}
                className={all ? 'disabled' : ''}
                block
              >
                Selected Aetiology
              </Button>
            </Col>
            <Col lg="2" md="3" sm="4" className="mb-3 d-flex">
              <Button
                size="sm"
                variant="primary"
                onClick={() => mergeAetiology({ all: true })}
                className={!all ? 'disabled' : ''}
                block
              >
                All Aetiologies
              </Button>
            </Col>
          </Row>
        </>
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
    if (data.length) {
      const aetiologies = [...new Set(data.map((obj) => obj.Aetiology))]
        .sort()
        .map((Aetiology) => (
          <Col key={Aetiology} lg="2" md="3" sm="4" className="mb-3 d-flex">
            <Button
              size="sm"
              variant="dark"
              onClick={() =>
                mergeAetiology({
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
        ));

      return (
        <>
          <Row className="justify-content-center">{aetiologies}</Row>
          <Row className="justify-content-center">
            <Col lg="2" md="3" sm="4" className="mb-3 d-flex">
              <Button
                size="sm"
                variant="primary"
                onClick={() => mergeAetiology({ all: false })}
                className={all ? 'disabled' : ''}
                block
              >
                Selected Aetiology
              </Button>
            </Col>
            <Col lg="2" md="3" sm="4" className="mb-3 d-flex">
              <Button
                size="sm"
                variant="primary"
                onClick={() => mergeAetiology({ all: true })}
                className={!all ? 'disabled' : ''}
                block
              >
                All Aetiologies
              </Button>
            </Col>
          </Row>
        </>
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
      thumbnails.forEach((t) => URL.revokeObjectURL(t));

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
                path: `msigportal/Database/Aetiology/Profile_logo/${Signature}.svg`,
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
      mergeAetiology({ thumbnails: newThumbnails });
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
              path: `msigportal/Database/Aetiology/Profile_logo/${tissue['Tissue Specific Signature']}.svg`,
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
      mergeAetiology({ tissueThumbnails: newThumbnails });
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
              path: `msigportal/Database/Aetiology/Profile_logo/${refSig['Ref Signature']}.svg`,
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
      mergeAetiology({ refSigThumbnails: newThumbnails });
    } catch (err) {
      mergeError(err.message);
    }
  }

  function standardView() {
    function getAllSignatures() {
      if (data.length) {
        return thumbnails.map(
          ({ Aetiology, Signature, url, Source_URL }, index) => (
            <Col key={index} lg="1" md="3" sm="4" className="mb-2 px-1">
              <div
                onClick={() =>
                  mergeAetiology({
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
      if (thumbnails.length) {
        return thumbnails
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
                    mergeAetiology({
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
      if (data.length) {
        return [
          ...new Set(
            data
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
              onClick={() => mergeAetiology({ study: Study })}
              className={study != Study ? 'disabled' : ''}
              block
            >
              {Study}
            </Button>
          </Col>
        ));
      }
    }

    function getInfo() {
      if (data.length && selectedSource) {
        const info = data.filter(
          ({ Source_URL, URL: alias_URL }) =>
            Source_URL == selectedSource || alias_URL == selectedSource
        )[0];

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
                info.Description.map((text) => <p>{text}</p>)
              ))}

            {profileURL ? (
              <div>
                <Plot
                  className="p-3 border rounded mb-3"
                  maxHeight={'300px'}
                  plotURL={profileURL}
                />
                <Row className="justify-content-center">{getStudy()}</Row>

                {exposureURL ? (
                  <>
                    <p>
                      Select the cancer study to review the TMB of selected
                      signatures. TMB shown as the numbers of mutations per
                      megabase (log10) attributed to each mutational signature
                      in samples where the signature is present. Only those
                      cancer types with tumors in which signature activity is
                      attributed are shown. The numbers below the dots for each
                      cancer type indicate the number of tumors in which the
                      signatures was attributed (above the horizontal bar, in
                      blue) and the total number of tumors analyzed (below the
                      blue bar, in green).
                    </p>
                    <Plot
                      className="p-3 border"
                      maxHeight={'400px'}
                      plotURL={exposureURL}
                    />
                  </>
                ) : (
                  <div className="p-3 border">
                    <p>
                      A signature was not detected in any sample of the selected
                      study
                    </p>
                  </div>
                )}
              </div>
            ) : (
              <div className="my-5">
                <LoadingOverlay active={true} />
              </div>
            )}
          </div>
        );
      }
    }

    return (
      <div>
        <div className="mx-auto p-3 pt-0">
          {getAetiologies()}
          <hr className="mb-3" />
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
        <div className="mx-auto px-5 py-3">{getInfo()}</div>
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
                  mergeAetiology({
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
                    mergeAetiology({
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
      if (refSigThumbnails.length) {
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
                  onClick={() => mergeAetiology({ refSig: v['Ref Signature'] })}
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
        return (
          <div className="my-5">
            <LoadingOverlay active={true} />
          </div>
        );
      }
    }

    function getCancerSpecificInfo() {
      if (data.length && tissue && refSig) {
        let info = data.filter(
          (v) =>
            v['Tissue Specific Signature'] == tissue &&
            v['Ref Signature'] == refSig
        )[0];

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
              info.Description.map((text) => <p>{text}</p>)
            )}

            <Plot
              className="p-3 border rounded mb-3"
              maxHeight={'300px'}
              plotURL={refSigURL}
            />

            <Plot
              className="p-3 border rounded mb-3"
              maxHeight={'300px'}
              plotURL={tissueURL}
            />

            <Plot
              className="p-3 border rounded mb-3"
              maxHeight={'300px'}
              plotURL={exposureURL}
            />
          </div>
        );
      }
    }
    return (
      <div>
        <div className="mx-auto p-3">
          {getCancerAetiology()}
          <hr className="mb-3" />
          <Row className={`justify-content-center ${all ? 'd-none' : ''}`}>
            {getTissues()}
          </Row>
          <Row
            className={`px-2 justify-content-center ${!all ? 'd-none' : ''}`}
          >
            {getAllTissues()}
          </Row>
          <hr className="mb-3" />
          <Row className="justify-content-center">{getRefSig()}</Row>
        </div>
        <hr />
        <div className="mx-auto px-5 py-3">{getCancerSpecificInfo()}</div>
      </div>
    );
  }

  return (
    <div className="bg-white border rounded">
      <div className="mx-auto pt-3">
        <Row className="justify-content-center">{getCategories()}</Row>
        <hr />
      </div>
      {category != 'Cancer Specific Signature' && standardView()}
      {category == 'Cancer Specific Signature' && cancerSpecificView()}
    </div>
  );
}
