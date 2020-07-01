import React, { useState, useEffect } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import {
  faChevronRight,
  faChevronLeft,
} from '@fortawesome/free-solid-svg-icons';
import { useSelector } from 'react-redux';
import {
  store,
  updateVisualizeResults,
  updateError,
} from '../../../services/store';

const { Group, Label, Control } = Form;

const root =
  process.env.NODE_ENV === 'development'
    ? 'http://localhost:8330/'
    : window.location.pathname;

export default function Results() {
  const {
    projectID,
    mapping,
    plots,
    displayedPlot,
    plotURL,
    error,
  } = useSelector((state) => state.visualizeResults);

  useEffect(() => {
    if (projectID.length) {
      getPlotMapping(projectID);
    }
  }, [projectID]);

  useEffect(() => {
    if (plots.length && !displayedPlot.length) {
      setPlot(mapping[0].Sample_Name);
    }
  }, [mapping]);

  async function getPlotMapping() {
    const response = await fetch(`${root}visualizeSummary`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ projectID: projectID }),
    });
    const mapping = await response.json();
    store.dispatch(
      updateVisualizeResults({
        mapping: mapping,
        plots: mapping.map((plot) => plot.Sample_Name).sort(),
      })
    );
  }

  async function setPlot(plotName) {
    const plotPath = mapping.filter((plot) => {
      return plot.Sample_Name == plotName;
    });
    const response = await fetch(`${root}svg`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ path: plotPath[0].Location }),
    });
    if (!response.ok) {
      const { msg } = await response.json();
      store.dispatch(updateError({ visible: true, message: msg }));
    } else {
      const pic = await response.blob();
      let objectURL = URL.createObjectURL(pic);
      store.dispatch(
        updateVisualizeResults({
          displayedPlot: plotName,
          plotURL: objectURL,
        })
      );
    }
  }

  function nextPlot() {
    const currIndex = plots.indexOf(displayedPlot);
    currIndex < plots.length - 1
      ? setPlot(plots[currIndex + 1])
      : setPlot(plots[0]);
  }

  function prevPlot() {
    const currIndex = plots.indexOf(displayedPlot);
    currIndex > 0
      ? setPlot(plots[currIndex - 1])
      : setPlot(plots[plots.length - 1]);
  }

  return error.length ? (
    <h4 className="text-danger">{error}</h4>
  ) : plots.length ? (
    <div>
      <Form>
        <Group>
          <Label>View Plots</Label>
          <Row className="justify-content-center">
            <Col sm="auto">
              <button className="faButton" onClick={() => prevPlot()}>
                <FontAwesomeIcon icon={faChevronLeft} size="2x" />
              </button>
            </Col>
            <Col sm="auto">
              <Control
                as="select"
                value={displayedPlot}
                onChange={(e) => setPlot(e.target.value)}
              >
                {plots.map((plot) => {
                  return (
                    <option key={plot} value={plot}>
                      {plot}
                    </option>
                  );
                })}
              </Control>
            </Col>
            <Col sm="auto">
              <button className="faButton" onClick={() => nextPlot()}>
                <FontAwesomeIcon icon={faChevronRight} size="2x" />
              </button>
            </Col>
          </Row>
        </Group>
      </Form>
      <img className="w-100" src={plotURL}></img>
    </div>
  ) : (
    <div>
      <h2>Instructions</h2>
      <p>Upload Sample Variants and specify parameters</p>
    </div>
  );
}
