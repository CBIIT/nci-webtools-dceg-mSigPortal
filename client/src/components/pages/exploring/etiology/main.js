import React, { useState } from 'react';
import { Form, Row, Col, Accordion, Card, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { dispatchError, dispatchExploring } from '../../../../services/store';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;
const { Group, Label, Check, Control } = Form;

export default function EtiologyExploring() {
  const rootURL = window.location.pathname;
  const { displayTab, exposureAccordion, publicDataOptions } = useSelector(
    (state) => state.exploring
  );
  const {
    study,
    studyOptions,
    strategy,
    strategyOptions,
    cancer,
    cancerOptions,
    refSigData,
    refSignatureSet,
    refSignatureSetOptions,
    genome,
    genomeOptions,
    exposureFile,
    matrixFile,
    signatureFile,
    usePublicSignature,
    source,
    loading,
    projectID,
  } = useSelector((state) => state.expExposure);

  const sections = [
    // {
    //   component: <Tumor />,
    //   id: 'tumor',
    //   title: 'Tumor Mutational Burden',
    // },
    // {
    //   component: <Activity calculateActivity={calculateActivity} />,
    //   id: 'activity',
    //   title: 'Mutational Signature Activity',
    // },
    // {
    //   component: <Decomposition />,
    //   id: 'decomposition',
    //   title: 'Evaluating the Performance of Mutational Signature Decomposition',
    // },
    // {
    //   component: <Association calculateAssociation={calculateAssociation} />,
    //   id: 'association',
    //   title: 'Mutational Signature Association',
    // },
    // {
    //   component: <Landscape calculateLandscape={calculateLandscape} />,
    //   id: 'landscape',
    //   title: 'Landscape of Mutational Signature Activity',
    // },
    // {
    //   component: <Prevalence calculatePrevalence={calculatePrevalence} />,
    //   id: 'prevalence',
    //   title: 'Prevalence of Mutational Signature',
    // },
  ];

  return (
    <div className="position-relative">
      <h4>TBA</h4>
    </div>
  );
}
