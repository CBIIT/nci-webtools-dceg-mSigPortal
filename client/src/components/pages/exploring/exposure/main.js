import React from 'react';
import { Form, Row, Col, Accordion, Card, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import Tumor from './tumor';
import Activity from './activity';
import Association from './association';
import Decomposition from './decomposition';
import Landscape from './landscape';
import Prevalence from './prevalence';
import {
  dispatchError,
  dispatchExploring,
  dispatchExpExposure,
  dispatchExpTumor,
  dispatchExpActivity,
  dispatchExpAssociation,
  dispatchExpDecomposition,
  dispatchExpLandscape,
  dispatchExpPrevalence,
} from '../../../../services/store';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function ExposureExploring() {
  const rootURL = window.location.pathname;
  const {
    displayTab,
    projectID,
    exposureAccordion,
    publicDataOptions,
  } = useSelector((state) => state.exploring);
  const {
    study,
    studyOptions,
    strategy,
    strategyOptions,
    refSigData,
    refSignatureSet,
    refSignatureSetOptions,
    datasource,
    loading,
  } = useSelector((state) => state.expExposure);
  const activityArgs = useSelector((state) => state.expActivity);
  const associationArgs = useSelector((state) => state.expAssociation);
  const landscapeArgs = useSelector((state) => state.expLandscape);
  const prevalenceArgs = useSelector((state) => state.expPrevalence);

  function submitR(fn, args) {
    return fetch(`${rootURL}exploringR`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        fn: fn,
        args: args,
        projectID: projectID,
      }),
    }).then((res) => res.json());
  }

  async function calculateActivity() {
    dispatchExpActivity({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const { debugR, output } = await submitR('exposurePublic', {
        fn: 'activity',
        common: JSON.stringify({
          study: study,
          strategy: strategy,
          refSignatureSet: refSignatureSet,
        }),
        activity: JSON.stringify({ signatureName: activityArgs.signatureName }),
      });

      if (output) {
        if (output.activityPath)
          dispatchExpActivity({
            plotPath: output.activityPath,
            debugR: debugR,
            err: false,
          });
        else dispatchExpActivity({ err: true, debugR: debugR });
      }
    } catch (err) {
      dispatchError(err);
    }

    dispatchExpActivity({ loading: false });
  }

  async function calculateAssociation() {
    dispatchExpAssociation({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const { debugR, output } = await submitR('exposurePublic', {
        fn: 'association',
        common: JSON.stringify({
          study: study,
          strategy: strategy,
          refSignatureSet: refSignatureSet,
        }),
        association: JSON.stringify({
          cancerType:
            associationArgs.cancer == 'None' ? null : associationArgs.cancer,
          both: associationArgs.both,
          signatureName1: associationArgs.signatureName1,
          signatureName2: associationArgs.signatureName2,
        }),
      });

      if (output) {
        if (output.associationPath)
          dispatchExpAssociation({
            plotPath: output.associationPath,
            debugR: debugR,
            err: false,
          });
        else dispatchExpAssociation({ err: true, debugR: debugR });
      }
    } catch (err) {
      dispatchError(err);
    }

    dispatchExpAssociation({ loading: false });
  }

  async function calculateLandscape() {
    dispatchExpLandscape({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const { debugR, output } = await submitR('exposurePublic', {
        fn: 'landscape',
        common: JSON.stringify({
          study: study,
          strategy: strategy,
          refSignatureSet: refSignatureSet,
        }),
        landscape: JSON.stringify({
          cancerType: landscapeArgs.cancer,
          varDataPath: landscapeArgs.varDataPath,
        }),
      });

      if (output) {
        if (output.landscapePath)
          dispatchExpLandscape({
            plotPath: output.landscapePath,
            debugR: debugR,
            err: false,
          });
        else dispatchExpLandscape({ err: true, debugR: debugR });
      }
    } catch (err) {
      dispatchError(err);
    }

    dispatchExpLandscape({ loading: false });
  }

  async function calculatePrevalence() {
    dispatchExpPrevalence({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const { debugR, output } = await submitR('exposurePublic', {
        fn: 'prevalence',
        common: JSON.stringify({
          study: study,
          strategy: strategy,
          refSignatureSet: refSignatureSet,
        }),
        prevalence: JSON.stringify({
          cancerType: prevalenceArgs.cancer,
          mutation: parseFloat(prevalenceArgs.mutation) || 100,
        }),
      });

      if (output) {
        if (output.prevalencePath)
          dispatchExpPrevalence({
            plotPath: output.prevalencePath,
            debugR: debugR,
            err: false,
          });
        else dispatchExpPrevalence({ err: true, debugR: debugR });
      }
    } catch (err) {
      dispatchError(err);
    }

    dispatchExpPrevalence({ loading: false });
  }

  async function handleCalculate(fn) {
    dispatchExpExposure({ loading: true });

    const { debugR, output } = await submitR('exposurePublic', {
      fn: fn,
      common: JSON.stringify({
        study: study,
        strategy: strategy,
        refSignatureSet: refSignatureSet,
      }),
      activity: JSON.stringify({ signatureName: activityArgs.signatureName }),
      association: JSON.stringify({
        cancerType:
          associationArgs.cancer == 'None' ? null : associationArgs.cancer,
        both: associationArgs.both,
        signatureName1: associationArgs.signatureName1,
        signatureName2: associationArgs.signatureName2,
      }),
      landscape: JSON.stringify({
        cancerType: landscapeArgs.cancer,
        varDataPath: landscapeArgs.varDataPath,
      }),
      prevalence: JSON.stringify({
        cancerType: prevalenceArgs.cancer,
        mutation: parseFloat(prevalenceArgs.mutation) || 100,
      }),
    });

    if (output) {
      if (output.tumorPath)
        dispatchExpTumor({
          plotPath: output.tumorPath,
          debugR: debugR,
          err: false,
        });
      else dispatchExpTumor({ err: true, debugR: debugR });
      if (output.activityPath)
        dispatchExpActivity({
          plotPath: output.activityPath,
          debugR: debugR,
          err: false,
        });
      else dispatchExpActivity({ err: true, debugR: debugR });
      if (output.associationPath)
        dispatchExpAssociation({
          plotPath: output.associationPath,
          debugR: debugR,
          err: false,
        });
      else dispatchExpAssociation({ err: true, debugR: debugR });
      if (output.decompositionPath)
        dispatchExpDecomposition({
          plotPath: output.decompositionPath,
          txtPath: output.decompositionData,
          debugR: debugR,
          err: false,
        });
      else dispatchExpDecomposition({ err: true, debugR: debugR });
      if (output.landscapePath)
        dispatchExpLandscape({
          plotPath: output.landscapePath,
          debugR: debugR,
          err: false,
        });
      else dispatchExpLandscape({ err: true, debugR: debugR });
      if (output.prevalencePath)
        dispatchExpPrevalence({
          plotPath: output.prevalencePath,
          debugR: debugR,
          err: false,
        });
      else dispatchExpPrevalence({ err: true, debugR: debugR });
    }

    dispatchExpExposure({ loading: false });
  }

  function handleStudy(study) {
    const strategyOptions = [
      ...new Set(
        publicDataOptions
          .filter((data) => data.Study == study)
          .map((data) => data.Dataset)
      ),
    ];

    const cancerOptions = [
      ...new Set(
        publicDataOptions
          .filter((data) => data.Study == study)
          .map((data) => data.Cancer_Type)
      ),
    ];

    dispatchExpExposure({
      study: study,
      strategy: strategyOptions[0],
      strategyOptions: strategyOptions,
    });

    dispatchExpLandscape({
      cancer: cancerOptions[0],
      cancerOptions: cancerOptions,
    });
  }

  function handleSet(set) {
    const signatureNameOptions = [
      ...new Set(
        refSigData
          .filter((row) => row.Signature_set_name == set)
          .map((row) => row.Signature_name)
      ),
    ];
    dispatchExpExposure({
      refSignatureSet: set,
      signatureNameOptions: signatureNameOptions,
    });
  }

  const sections = [
    {
      component: <Tumor />,
      id: 'tumor',
      title: 'Tumor Mutational Burden',
    },
    {
      component: <Activity calculateActivity={calculateActivity} />,
      id: 'activity',
      title: 'Mutational Signature Activity',
    },
    {
      component: <Association calculateAssociation={calculateAssociation} />,
      id: 'association',
      title: 'Mutational Signature Association',
    },
    {
      component: <Decomposition />,
      id: 'decomposition',
      title: 'Evaluating the Performance of Mutational Signature Decomposition',
    },
    {
      component: <Landscape calculateLandscape={calculateLandscape} />,
      id: 'landscape',
      title: 'Landscape of Mutational Signature Activity',
    },
    {
      component: <Prevalence calculatePrevalence={calculatePrevalence} />,
      id: 'prevalence',
      title: 'Prevalence of Mutational Signature',
    },
  ];
  return (
    <div className="position-relative">
      <Form>
        <LoadingOverlay active={loading} />
        <div>
          <Row className="justify-content-center">
            <Col sm="2">
              <Select
                id="tumorStudy"
                label="Study"
                value={study}
                options={studyOptions}
                onChange={handleStudy}
              />
            </Col>
            <Col sm="2">
              <Select
                id="tumorStrategy"
                label="Experimental Strategy"
                value={strategy}
                options={strategyOptions}
                onChange={(strategy) =>
                  dispatchExpExposure({ strategy: strategy })
                }
              />
            </Col>
            <Col sm="4">
              <Select
                id="tumorSet"
                label="Reference Signature Set"
                value={refSignatureSet}
                options={refSignatureSetOptions}
                onChange={handleSet}
              />
            </Col>
            <Col sm="3"></Col>
            <Col sm="1" className="m-auto">
              <Button variant="primary" onClick={() => handleCalculate('all')}>
                Calculate
              </Button>
            </Col>
          </Row>
        </div>
      </Form>
      {sections.map(({ component, id, title }) => {
        return (
          <Accordion activeKey={exposureAccordion[id]} key={id}>
            <Card>
              <Toggle
                className="font-weight-bold"
                as={Header}
                eventKey={exposureAccordion[id]}
                onClick={() =>
                  dispatchExploring({
                    exposureAccordion: {
                      ...exposureAccordion,
                      [id]: !exposureAccordion[id],
                    },
                  })
                }
              >
                {exposureAccordion[id] == true ? (
                  <FontAwesomeIcon icon={faMinus} />
                ) : (
                  <FontAwesomeIcon icon={faPlus} />
                )}{' '}
                {title}
              </Toggle>
              <Collapse eventKey={exposureAccordion[id]}>
                <Body>{component}</Body>
              </Collapse>
            </Card>
          </Accordion>
        );
      })}
    </div>
  );
}
