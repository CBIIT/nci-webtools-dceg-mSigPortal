import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
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
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;
const root =
  process.env.NODE_ENV === 'development'
    ? 'http://localhost:8330/'
    : window.location.pathname;

export default function Results() {
  const [downloadOverlay, setOverlay] = useState(false);
  const {
    projectID,
    mapping,
    displayedPlotIndex,
    plotURL,
    error,
    debug,
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
    const response = await fetch(`${root}visualize/summary`, {
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
      const response = await fetch(`${root}visualize/svg`, {
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

  function getPlotName() {
    if (displayedPlotIndex) {
      const plot = mapping[displayedPlotIndex];
      return `${plot.Sample_Name}-${plot.Profile_Type}-${plot.Matrix}.svg`;
    } else {
      // assume index is 0 upon loading a projectID
      const plot = mapping[0];
      return `${plot.Sample_Name}-${plot.Profile_Type}-${plot.Matrix}.svg`;
    }
  }

  async function downloadResults() {
    setOverlay(true);
    const response = await fetch(`${root}visualize/results`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ projectID: projectID }),
    });
    if (!response.ok) {
      setOverlay(false);
      const { msg } = await response.json();
      store.dispatch(updateError({ visible: true, message: msg }));
    } else {
      setOverlay(false);
      const req = await response.json();
      const id = req.projectID;
      const tempLink = document.createElement('a');

      tempLink.href = `${root}visualize/download?id=${id}`;
      document.body.appendChild(tempLink);
      tempLink.click();
      document.body.removeChild(tempLink);
    }
  }

  return error.length ? (
    <h4 className="text-danger">{error}</h4>
  ) : mapping.length ? (
    <div>
      <Form>
        <Group controlId="selectPlot">
          <span className="d-flex">
            <Label className="px-2 py-1">Results</Label>
          </span>
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
      <Row>
        <Col>
          <img className="w-100 my-4" src={plotURL}></img>
        </Col>
      </Row>
      <div className="d-flex">
        <a className="ml-2 px-2 py-1" href={plotURL} download={getPlotName()}>
          Download Plot
        </a>
        <span className="ml-auto">
          <LoadingOverlay
            active={downloadOverlay}
            // showIndicator={true}
          />
          <Button
            className="px-2 py-1"
            variant="link"
            onClick={() => downloadResults()}
          >
            Download Results
          </Button>
        </span>
      </div>
      <div>
        <div>stdout</div>
        <pre>{debug.stdout}</pre>
        <div>stderr</div>
        <pre>{debug.stderr}</pre>
      </div>
    </div>
  ) : (
    <div>
      <h2>Instructions</h2>
      <p>Upload Sample Variants and specify parameters</p>
    </div>
  );
}
