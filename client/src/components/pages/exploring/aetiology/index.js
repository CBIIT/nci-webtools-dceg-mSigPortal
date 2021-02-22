import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Accordion, Card, Button, Nav } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';

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
            await fetch(`api/public/Others/json/aetiology.json`)
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
      return [...new Set(data.map((obj) => obj.aetiology))];
    } else {
      return [];
    }
  }

  function getSignatures() {
    if (data.length) {
      return [...new Set(data.map((obj) => obj.signature))];
    } else {
      return [];
    }
  }

  return (
    <div className="bg-white border rounded p-4">
      <div className="mx-auto" style={{ width: '300px' }}>
        <Row className="d-flex justify-content-center">
          {getAetiologies().map((aetiology) => (
            <Col md="4" className="mb-3">
              <Button
                className="d-flex mx-auto"
                onClick={() => mergeAetiology({ aetiology: aetiology })}
              >
                {aetiology}
              </Button>
            </Col>
          ))}
        </Row>
        <Row className="d-flex justify-content-center">
          {getSignatures().map((signature) => (
            <Col md="4" className="mb-3">
              <Button
                className="d-flex mx-auto"
                onClick={() => mergeAetiology({ signature: signature })}
              >
                {signature}
              </Button>
            </Col>
          ))}
        </Row>
      </div>
      {/* <pre>
        <code>
          {JSON.stringify(
            (() => {
              const { data, ...rest } = exploring.aetiology;
              return rest;
            })()
          )}
        </code>
      </pre> */}
    </div>
  );
}
