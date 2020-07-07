import React, { useEffect } from 'react';
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
    displayedPlotIndex,
    plotURL,
    error,
  } = useSelector((state) => state.visualizeResults);

  useEffect(() => {
    if (projectID.length) {
      getPlotMapping(projectID);
    }
  }, [projectID]);

  // load first plot after results are recieved
  useEffect(() => {
    if (mapping.length && !displayedPlotIndex.length) {
      setPlot(0);
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
      })
    );
  }

  async function setPlot(index) {
    const plot = mapping[index];
    if (plot) {
      const response = await fetch(`${root}svg`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ path: plot.Location }),
      });
      if (!response.ok) {
        const { msg } = await response.json();
        store.dispatch(updateError({ visible: true, message: msg }));
      } else {
        const pic = await response.blob();
        const objectURL = URL.createObjectURL(pic);

        if (plotURL.length) URL.revokeObjectURL(plotURL);
        store.dispatch(
          updateVisualizeResults({
            displayedPlotIndex: index,
            plotURL: objectURL,
          })
        );
      }
    } else {
      console.log('invalid index', index);
    }
  }

  function nextPlot() {
    displayedPlotIndex < mapping.length - 1
      ? setPlot(parseInt(displayedPlotIndex + 1))
      : setPlot(0);
  }

  function prevPlot() {
    displayedPlotIndex > 0
      ? setPlot(parseInt(displayedPlotIndex - 1))
      : setPlot(mapping.length - 1);
  }

  function getDownloadName() {
    if (displayedPlotIndex) {
      const plot = mapping[displayedPlotIndex];
      return `${plot.Sample_Name}-${plot.Profile_Type}-${plot.Matrix}.svg`;
    } else {
      // assume index is 0 upon loading a projectID
      const plot = mapping[0];
      return `${plot.Sample_Name}-${plot.Profile_Type}-${plot.Matrix}.svg`;
    }
  }

  return error.length ? (
    <h4 className="text-danger">{error}</h4>
  ) : mapping.length ? (
    <div>
      <Form>
        <Group controlId="selectPlot">
          <Label>View Plots</Label>
          <Row className="justify-content-center">
            <Col sm="auto">
              <button
                className="faButton navButton"
                onClick={(e) => {
                  e.preventDefault();
                  prevPlot();
                }}
              >
                <FontAwesomeIcon icon={faChevronLeft} size="2x" />
              </button>
            </Col>
            <Col sm="auto">
              <Control
                as="select"
                value={displayedPlotIndex}
                onChange={(e) => setPlot(e.target.value)}
              >
                {mapping.map((plot, index) => {
                  return (
                    <option key={index} value={index}>
                      {`${plot.Sample_Name} | ${plot.Profile_Type} | ${plot.Matrix}`}
                    </option>
                  );
                })}
              </Control>
            </Col>
            <Col sm="auto">
              <button
                className="faButton navButton"
                onClick={(e) => {
                  e.preventDefault();
                  nextPlot();
                }}
              >
                <FontAwesomeIcon icon={faChevronRight} size="2x" />
              </button>
            </Col>
          </Row>
        </Group>
      </Form>
      <a href={plotURL} download={getDownloadName()}>
        Download Plot
      </a>
      <img className="w-100" src={plotURL}></img>
    </div>
  ) : (
    <div>
      <h2>Instructions</h2>
      <p>Upload Sample Variants and specify parameters</p>
    </div>
  );
}
