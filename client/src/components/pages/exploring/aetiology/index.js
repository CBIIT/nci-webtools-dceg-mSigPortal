import React, { useState, useEffect } from 'react';
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
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeAetiology = (state) =>
    dispatch(actions.mergeExploring({ aetiology: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    category,
    aetiology,
    signature,
    study,
    all,
    data,
    thumbnails,
    sigURL,
    tmbURL,
  } = exploring.aetiology;

  const categories = [
    { name: 'Cosmic Mutational Signatures', file: 'Aetiology_cosmic.json' },
    {
      name: 'Environmental Mutagenesis',
      file: 'Aetiology_enviromental_mutagenesis.json',
    },
    { name: 'Gene Edits', file: 'Aetiology_gene_edits.json' },
    { name: 'Cancer Specific Signature', file: '' },
    { name: 'Others', file: '' },
  ];

  // useEffect(() => {
  //   mergeExploring({ displayTab: 'aetiology' });
  // }, []);

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
      if (sigURL) URL.revokeObjectURL(sigURL);
      if (tmbURL) URL.revokeObjectURL(tmbURL);
      const [sig, tmb] = await Promise.all([
        getImageS3(`msigportal/Database/Aetiology/Profile/${signature}.svg`),
        getImageS3(
          `msigportal/Database/Aetiology/Exposure/${signature}_${study}.svg`
        ),
      ]);
      mergeAetiology({ sigURL: sig, tmbURL: tmb });
    };

    if (signature && study) getPlots();
  }, [signature, study]);

  // get thumbnails
  useEffect(() => {
    if (data.length) {
      const signatures = data
        .slice()
        .filter(({ Study }) => Study == study)
        .sort(naturalSort)
        .sort(profileSort);

      getThumbnails(signatures);
    }
  }, [data]);

  function getCategories() {
    return categories.map(({ name, file }) => (
      <Col key={name} lg="2" md="3" sm="4" className="mb-3 d-flex">
        <Button
          size="sm"
          variant="dark"
          className="d-flex mx-auto"
          onClick={
            file && name != category
              ? async () => {
                  mergeAetiology({
                    category: name,
                    aetiology: '',
                    signature: '',
                    study: '',
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
      const aetiologies = [...new Set(data.map((obj) => obj.Aetiology))].map(
        (Aetiology) => (
          <Col key={Aetiology} lg="2" md="3" sm="4" className="mb-3 d-flex">
            <Button
              size="sm"
              variant="dark"
              className="d-flex mx-auto"
              onClick={() =>
                mergeAetiology({
                  aetiology: Aetiology,
                  signature: '',
                })
              }
              className={aetiology != Aetiology ? 'disabled' : ''}
              block
            >
              {Aetiology}
            </Button>
          </Col>
        )
      );

      return (
        <>
          <Row className="justify-content-center">{aetiologies}</Row>
          <Row className="justify-content-center">
            <Col lg="2" md="3" sm="4" className="mb-3 d-flex">
              <Button
                size="sm"
                variant="primary"
                className="d-flex mx-auto"
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
                className="d-flex mx-auto"
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
      return [];
    }
  }

  function profileSort(a, b) {
    const sigOrder = ['SBS', 'DBS', 'ID'];
    let c = 0,
      d = 0;

    sigOrder.forEach((profile, i) => {
      if (a.Signature.includes(profile)) {
        c = i;
      }
      if (b.Signature.includes(profile)) {
        d = i;
      }
    });

    return c - d;
  }

  function naturalSort(a, b) {
    return a.Signature.localeCompare(b.Signature, undefined, {
      numeric: true,
      sensitivity: 'base',
    });
  }

  async function getThumbnails(signatures) {
    try {
      thumbnails.forEach((t) => URL.revokeObjectURL(t));

      const newThumbnails = await Promise.all(
        signatures.map(({ Aetiology, Study, Signature }) =>
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

  function getAllSignatures() {
    if (data.length) {
      return thumbnails.map(({ Aetiology, Signature, url }, i) => (
        <Col key={Signature + i} lg="1" md="3" sm="4" className="mb-2 px-1">
          <div
            onClick={() =>
              mergeAetiology({
                aetiology: Aetiology,
                signature: Signature,
              })
            }
            className={`sigIcon border rounded ${
              aetiology != Aetiology
                ? 'inactive'
                : signature == Signature
                ? 'active'
                : ''
            }`}
            title={`${Aetiology} - ${Signature}`}
          >
            <img
              src={url}
              className="w-100"
              // height="70"
              alt=""
              // alt={Signature}
            />
            <strong className="sigLabel">{Signature}</strong>
          </div>
        </Col>
      ));
    } else {
      return [];
    }
  }

  function getSignatures() {
    if (thumbnails.length) {
      return thumbnails
        .filter(({ Aetiology }) => Aetiology == aetiology)
        .map(({ Signature, url }) => {
          return (
            <Col key={Signature} md="2" sm="4" className="mb-3">
              <div
                className={`sigIcon border rounded ${
                  signature == Signature ? 'active' : ''
                }`}
                title={`${aetiology} - ${Signature}`}
                onClick={() => mergeAetiology({ signature: Signature })}
              >
                <img
                  src={url}
                  className="w-100"
                  // height="110"
                  alt=""
                  // alt={Signature}
                />
                <strong className="sigLabel">{Signature}</strong>
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
        <Col lg="2" md="3" sm="4" className="mb-3 d-flex">
          <Button
            size="sm"
            variant="primary"
            className="d-flex mx-auto"
            onClick={() => mergeAetiology({ study: Study })}
            className={study != Study ? 'disabled' : ''}
            block
          >
            {Study}
          </Button>
        </Col>
      ));
    } else {
      return (
        <div className="my-5">
          <LoadingOverlay active={true} />
        </div>
      );
    }
  }

  function getInfo() {
    if (data.length && signature) {
      const info = data.filter(
        ({ Study, Aetiology, Signature }) =>
          Study == study && Aetiology == aetiology && Signature == signature
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

          {sigURL ? (
            <div>
              <Plot
                className="p-3 border rounded mb-3"
                maxHeight={'300px'}
                plotURL={sigURL}
              />
              <Row className="justify-content-center">{getStudy()}</Row>
              <p>
                Select the cancer study to review the TMB of selected
                signatures. TMB shown as the numbers of mutations per megabase
                (log10) attributed to each mutational signature in samples where
                the signature is present. Only those cancer types with tumors in
                which signature activity is attributed are shown. The numbers
                below the dots for each cancer type indicate the number of
                tumors in which the signatures was attributed (above the
                horizontal bar, in blue) and the total number of tumors analyzed
                (below the blue bar, in green).
              </p>
              {tmbURL ? (
                <Plot
                  className="p-3 border"
                  maxHeight={'400px'}
                  plotURL={tmbURL}
                />
              ) : (
                <div className="p-3 border">
                  <p>
                    This signature was not detected in any sample of the
                    selected study
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
    <div className="bg-white border rounded">
      <div className="mx-auto p-3">
        <Row className="justify-content-center">{getCategories()}</Row>
        {getAetiologies()}
        <Row className={`justify-content-center ${all ? 'd-none' : ''}`}>
          {getSignatures()}
        </Row>
        <Row className={`px-2 justify-content-center ${!all ? 'd-none' : ''}`}>
          {getAllSignatures()}
        </Row>
      </div>
      <hr />
      <div className="mx-auto px-5 py-3">{getInfo()}</div>
    </div>
  );
}
