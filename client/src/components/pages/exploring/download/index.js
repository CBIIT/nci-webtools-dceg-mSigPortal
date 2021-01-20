import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Accordion, Card, Button, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { dispatchError, dispatchExploring } from '../../../../services/store';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;
const { Group, Label, Check, Control } = Form;

export default function Download() {
  useEffect(() => {
    dispatchExploring({ displayTab: 'download' });
  }, []);

  return <div className="bg-white border rounded p-4">TBA</div>;
}
