import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Card, Button } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import Select from '../../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
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

  const { aetiology, signature, study, all, data } = exploring.aetiology;

  useEffect(() => {
    mergeExploring({ displayTab: 'aetiology' });
  }, []);

  useEffect(() => {
    const getData = async () => {
      try {
        mergeAetiology({
          data: await (
            await fetch(`api/public/Aetiology/Aetiology.json`)
          ).json(),
        });
      } catch (_) {
        mergeError('Could not fetch Aetiology Info');
      }
    };

    if (!data.length) getData();
  }, [data]);

  // check if plot exists
  useEffect(() => {
    const checkPlot = async () => {
      let check = await fetch(
        `api/public/Aetiology/Exposure/${signature}_${study}.svg`,
        { method: 'HEAD' }
      );
      check.status === 200 ? setCheck(true) : setCheck(false);
    };

    if (signature && study) checkPlot();
  }, [signature, study]);

  const [checkPlot, setCheck] = useState(false);

  function getAetiologies() {
    if (data.length) {
      const aetiologies = [
        ...new Set(
          data.filter((obj) => obj.Study == study).map((obj) => obj.Aetiology)
        ),
      ].map((Aetiology) => (
        <Col lg="2" md="3" sm="4" className="mb-3 d-flex">
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
      ));

      return (
        <>
          {aetiologies}
          <Col lg="2" md="3" sm="4" className="mb-3 d-flex">
            <Button
              size="sm"
              variant="dark"
              className="d-flex mx-auto"
              onClick={() => mergeAetiology({ all: !all })}
              className={!all ? 'disabled' : ''}
              block
            >
              Toggle All Aetiologies
            </Button>
          </Col>
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

  function getAllSignatures() {
    if (data.length) {
      return data
        .filter(({ Study }) => Study == study)
        .sort(naturalSort)
        .sort(profileSort)
        .map(({ Aetiology, Signature }) => (
          <Col lg="1" md="3" sm="4" className="mb-3">
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
                src={`api/public/Aetiology/Profile_logo/${Signature}.svg`}
                className="w-100"
                // height="70"
                alt={Signature}
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
    if (data.length && aetiology != 'all') {
      const signatures = data
        .filter(
          ({ Study, Aetiology }) => Study == study && Aetiology == aetiology
        )
        .sort(naturalSort)
        .sort(profileSort);

      return signatures.map(({ Aetiology, Signature }) => (
        <Col md="2" sm="4" className="mb-3">
          <div
            className={`sigIcon border rounded ${
              aetiology != Aetiology
                ? 'inactive'
                : signature == Signature
                ? 'active'
                : ''
            }`}
            onClick={() => mergeAetiology({ signature: Signature })}
          >
            <img
              src={`api/public/Aetiology/Profile_logo/${Signature}.svg`}
              className="w-100"
              // height="110"
              alt={Signature}
            />
            <strong className="sigLabel">{Signature}</strong>
          </div>
        </Col>
      ));
    } else {
      return [];
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
            variant="dark"
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
      return [];
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
          {info.Source && (
            <div>
              <strong>Source: </strong>
              <a href={info.URL} target="_blank" rel="noreferrer">
                {info.Source}
              </a>
            </div>
          )}
          {typeof info.Description == 'string' ? (
            <p>{info.Description}</p>
          ) : (
            info.Description.map((text) => <p>{text}</p>)
          )}

          <Plot
            className="p-3 border rounded mb-3"
            maxHeight={'300px'}
            plotURL={`api/public/Aetiology/Profile/${signature}.svg`}
          />
          <Row className="justify-content-center">{getStudy()}</Row>
          {checkPlot ? (
            <Plot
              className="p-3 border"
              maxHeight={'400px'}
              plotURL={`api/public/Aetiology/Exposure/${signature}_${study}.svg`}
            />
          ) : (
            <div className="p-3 border">
              <p>
                This signature was not detected in any sample of the selected
                study
              </p>
            </div>
          )}
        </div>
      );
    }
  }

  return (
    <div className="bg-white border rounded">
      <div className="mx-auto p-3">
        <Row className="justify-content-center">{getAetiologies()}</Row>
        <Row className={`justify-content-center ${all ? 'd-none' : ''}`}>
          {getSignatures()}
        </Row>
        <Row className={`justify-content-center ${!all ? 'd-none' : ''}`}>
          {getAllSignatures()}
        </Row>
      </div>
      <hr />
      <div className="mx-auto px-5 py-3">{getInfo()}</div>
    </div>
  );
}
