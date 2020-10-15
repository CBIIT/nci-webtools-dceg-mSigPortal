import React from 'react';
import { Form, Row, Col, Accordion, Card, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import Tumor from './tumor';
import Activity from './activity';
import Decomposition from './decomposition';
import Landscape from './landscape';
import Prevalence from './prevalence';
import {
  dispatchError,
  dispatchExploring,
  dispatchExpExposure,
} from '../../../../services/store';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;
const { Group, Label, Control, Text } = Form;

export default function ExposureExploring() {
  const rootURL = window.location.pathname;
  const { displayTab, projectID, exposureAccordion } = useSelector(
    (state) => state.exploring
  );
  const {
    publicDataOptions,
    study,
    studyOptions,
    strategy,
    strategyOptions,
    refSignatureSet,
    refSignatureSetOptions,
    genomeSize,
    loading,
  } = useSelector((state) => state.expExposure);

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
    });
  }

  async function calculateTumor(fn, args) {
    console.log(fn);
    // dispatchExpTumor({
    //   loading: true,
    //   err: false,
    //   debugR: '',
    // });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        // dispatchExpTumor({
        //   loading: false,
        //   debugR: err,
        // });
      } else {
        const { debugR, output } = await response.json();

        // dispatchExpTumor({
        //   debugR: debugR,
        //   loading: false,
        //   plotPath: output.plotPath,
        //   txtPath: output.txtPath,
        // });
        // setRPlot(output.plotPath);
      }
    } catch (err) {
      dispatchError(err);
      // dispatchExpTumor({ loading: false });
    }
  }

  async function handleCalculate() {}

  function handleStudy(study) {
    const strategyOptions = [
      ...new Set(
        publicDataOptions
          .filter((data) => data.Study == study)
          .map((data) => data.Dataset)
      ),
    ];

    dispatchExpExposure({
      study: study,
      strategy: strategyOptions[0],
      strategyOptions: strategyOptions,
    });
  }

  const sections = [
    {
      component: <Tumor submitR={(fn, args) => submitR(fn, args)} />,
      id: 'tumor',
      title: 'Tumor Mutational Burden',
    },
    {
      component: <Activity submitR={(fn, args) => submitR(fn, args)} />,
      id: 'activity',
      title: 'Mutational Signature Activity',
    },
    {
      component: <Decomposition submitR={(fn, args) => submitR(fn, args)} />,
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
                onChange={(set) =>
                  dispatchExpExposure({ refSignatureSet: set })
                }
              />
            </Col>
            <Col sm="3">
              <Group controlId="tumorGenomeSize">
                <Label>Genome Size</Label>
                <Control
                  value={genomeSize}
                  onChange={(e) => {
                    dispatchExpExposure({
                      genomeSize: e.target.value,
                    });
                  }}
                />
                {/* <Text className="text-muted">(Ex. NCG>NTG)</Text> */}
              </Group>
            </Col>
            <Col sm="1" className="m-auto">
              <Button
                variant="primary"
                onClick={() => {
                  handleCalculate('tumorBurden', {
                    study: study,
                    strategy: strategy,
                    refSignatureSet: refSignatureSet,
                    genomeSize: parseFloat(genomeSize),
                  });
                }}
              >
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
