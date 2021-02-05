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
const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;
const { Group, Label, Check, Control } = Form;

export default function Download() {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeError = (state) => dispatch(actions.mergeModal({ error: state }));

  useEffect(() => {
    mergeExploring({ displayTab: 'download' });
  }, []);

  return <div className="bg-white border rounded p-4">TBA</div>;
}
