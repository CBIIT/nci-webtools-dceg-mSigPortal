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

  const { aetiology, signature, study, data } = exploring.aetiology;

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

  function getAetiologies() {
    if (data.length) {
      return [
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
              mergeAetiology({ aetiology: Aetiology, signature: '' })
            }
            className={aetiology != Aetiology ? 'disabled' : ''}
            block
          >
            {Aetiology}
          </Button>
        </Col>
      ));
    } else {
      return [];
    }
  }

  function getSignatures() {
    if (data.length) {
      return data
        .filter(({ Study }) => Study == study)
        .map(({ Aetiology, Signature }) => (
          <Col lg="2" md="3" sm="4" className="mb-3">
            <div
              className={`sigIcon border rounded ${
                aetiology != Aetiology
                  ? 'inactive'
                  : signature == Signature
                  ? 'active'
                  : ''
              }`}
            >
              <img
                src={`api/public/Aetiology/Profile_logo/${Signature}.svg`}
                className="w-100"
                onClick={() =>
                  mergeAetiology({
                    aetiology: Aetiology,
                    signature: Signature,
                  })
                }
                height="110"
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
      return (
        <select
          className="mb-3"
          onChange={(e) => mergeAetiology({ study: e.target.value })}
        >
          {[...new Set(data.map((obj) => obj.Study))].map((Study) => (
            <option value={Study}>{Study}</option>
          ))}
        </select>
      );
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
              <a href={info.SourceURL} target="_blank" rel="noreferrer">
                {info.Source}
              </a>
            </div>
          )}
          <p>{info.Description}</p>
          {/* {info.Description && info.Description.map((text) => <p>{text}</p>)} */}

          <Plot
            className="p-3 border"
            maxHeight={'500px'}
            plotURL={`api/public/Aetiology/Profile/${signature}.svg`}
          />
          <Plot
            className="p-3 border"
            maxHeight={'500px'}
            plotURL={`api/public/Aetiology/Exposure/${signature}_${study}.svg`}
          />
        </div>
      );
    }
  }

  return (
    <div className="bg-white border rounded">
      <div className="mx-auto p-3">
        <Row className="justify-content-center">{getStudy()}</Row>
        <Row className="justify-content-center">{getAetiologies()}</Row>
        <Row className="justify-content-center">{getSignatures()}</Row>
      </div>
      <hr />
      <div className="mx-auto px-5 py-3">{getInfo()}</div>
    </div>
  );
}
