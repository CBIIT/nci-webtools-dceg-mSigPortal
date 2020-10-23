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
} from '../../../../services/store';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;
const { Group, Label, Control, Text } = Form;

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

  // async function calculateTumor(fn, args) {
  //   console.log(fn);
  //   // dispatchExpTumor({
  //   //   loading: true,
  //   //   err: false,
  //   //   debugR: '',
  //   // });

  //   try {
  //     const response = await submitR(fn, args);
  //     if (!response.ok) {
  //       const err = await response.json();

  //       // dispatchExpTumor({
  //       //   loading: false,
  //       //   debugR: err,
  //       // });
  //     } else {
  //       const { debugR, output } = await response.json();

  //       // dispatchExpTumor({
  //       //   debugR: debugR,
  //       //   loading: false,
  //       //   plotPath: output.plotPath,
  //       //   txtPath: output.txtPath,
  //       // });
  //       // setRPlot(output.plotPath);
  //     }
  //   } catch (err) {
  //     dispatchError(err);
  //     // dispatchExpTumor({ loading: false });
  //   }
  // }

  async function handleCalculate() {
    dispatchExpExposure({ loading: true });

    const { debugR, output } = await submitR('exposurePublic', {
      common: JSON.stringify({
        study: study,
        strategy: strategy,
        refSignatureSet: refSignatureSet,
      }),
      activity: JSON.stringify({ signatureName: activityArgs.signatureName }),
      association: JSON.stringify({
        signatureName1: associationArgs.signatureName1,
        signatureName2: associationArgs.signatureName2,
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
      component: <Activity submitR={(fn, args) => submitR(fn, args)} />,
      id: 'activity',
      title: 'Mutational Signature Activity',
    },
    {
      component: <Association submitR={(fn, args) => submitR(fn, args)} />,
      id: 'association',
      title: 'Mutational Signature Association',
    },
    {
      component: <Decomposition />,
      id: 'decomposition',
      title: 'Evaluating the Performance of Mutational Signature Decomposition',
    },
    {
      component: <Landscape submitR={(fn, args) => submitR(fn, args)} />,
      id: 'landscape',
      title: 'Landscape of Mutational Signature Activity',
    },
    {
      component: <Prevalence submitR={(fn, args) => submitR(fn, args)} />,
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
              <Button variant="primary" onClick={() => handleCalculate()}>
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
