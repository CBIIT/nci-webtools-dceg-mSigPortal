import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button, Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualizeResults,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;

const root =
  process.env.NODE_ENV === 'development'
    ? 'http://localhost:8330/'
    : window.location.pathname;

export default function PyTab({ setPlot }) {
  const [downloadOverlay, setOverlay] = useState(false);
  const { projectID, pyPlotURL, pyTab } = useSelector(
    (state) => state.visualizeResults
  );

  const {
    mapping,
    filtered,
    selectName,
    selectProfile,
    selectMatrix,
    selectTag,
    nameOptions,
    profileOptions,
    matrixOptions,
    tagOptions,
    debugPy,
  } = pyTab;

  // set inital plot
  useEffect(() => {
    if (!pyPlotURL.length) {
      setPlot(0);
    }
  }, [filtered]);

  // change plot
  useEffect(() => {
    if (filtered.length && pyPlotURL.length) {
      setPlot(0);
    }
  }, [selectName, selectProfile, selectMatrix, selectTag]);

  function getPlotName() {
    const plot = filtered[0];
    return `${plot.Sample_Name}-${plot.Profile_Type}-${plot.Matrix}-${plot.Tag}.svg`;
  }

  async function downloadResults() {
    setOverlay(true);
    try {
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
        dispatchError(msg);
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
    } catch (err) {
      dispatchError(err);
    }
  }

  function filterSampleName(name) {
    const filteredPlots = mapping.filter((plot) => plot.Sample_Name == name);
    const profileOptions = [
      ...new Set(filteredPlots.map((plot) => plot.Profile_Type)),
    ];
    const matrixOptions = [
      ...new Set(
        filteredPlots
          .filter((plot) => plot.Profile_Type == profileOptions[0])
          .map((plot) => plot.Matrix)
      ),
    ];
    const tagOptions = [
      ...new Set(
        filteredPlots
          .filter((plot) => plot.Matrix == matrixOptions[0])
          .map((plot) => plot.Tag)
      ),
    ];

    dispatchVisualizeResults({
      pyTab: {
        ...pyTab,
        selectName: name,
        selectProfile: profileOptions[0],
        selectMatrix: matrixOptions[0],
        selectTag: tagOptions[0],
        profileOptions: profileOptions,
        matrixOptions: matrixOptions,
        tagOptions: tagOptions,
        filtered: filteredPlots,
      },
    });
  }

  function filterProfileType(profile) {
    const filteredPlots = mapping.filter(
      (plot) => plot.Profile_Type == profile
    );
    const matrixOptions = [
      ...new Set(filteredPlots.map((plot) => plot.Matrix)),
    ];
    const tagOptions = [
      ...new Set(
        filteredPlots
          .filter((plot) => plot.Matrix == matrixOptions[0])
          .map((plot) => plot.Tag)
      ),
    ];

    dispatchVisualizeResults({
      pyTab: {
        ...pyTab,
        selectProfile: profile,
        selectMatrix: matrixOptions[0],
        selectTag: tagOptions[0],
        matrixOptions: matrixOptions,
        tagOptions: tagOptions,
        filtered: filteredPlots,
      },
    });
  }

  function filterMatrix(matrix) {
    const filteredPlots = mapping.filter((plot) => plot.Matrix == matrix);
    const tagOptions = [...new Set(filteredPlots.map((plot) => plot.Tag))];

    dispatchVisualizeResults({
      pyTab: {
        ...pyTab,
        selectMatrix: matrix,
        selectTag: tagOptions[0],
        tagOptions: tagOptions,
        filtered: filteredPlots,
      },
    });
  }

  function filterTag(tag) {
    const filteredPlots = mapping.filter((plot) => plot.Tag == tag);

    dispatchVisualizeResults({
      pyTab: {
        ...pyTab,
        selectTag: tag,
        filtered: filteredPlots,
      },
    });
  }

  return (
    <div>
      <Form>
        <Group controlId="selectPlot">
          <Label className="py-1">Results</Label>
          <div className="border rounded p-2">
            <Row className="justify-content-center">
              <Col sm="3">
                <Label>Sample Name</Label>
                <Control
                  as="select"
                  value={selectName}
                  onChange={(e) => filterSampleName(e.target.value)}
                  custom
                >
                  <option value="0" disabled>
                    Select
                  </option>
                  {nameOptions.map((sampleName, index) => {
                    return (
                      <option key={index} value={sampleName}>
                        {sampleName}
                      </option>
                    );
                  })}
                </Control>
              </Col>
              <Col sm="3">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={selectProfile}
                  onChange={(e) => filterProfileType(e.target.value)}
                  custom
                >
                  <option value="0" disabled>
                    Select
                  </option>
                  {profileOptions.map((profile, index) => {
                    return (
                      <option key={index} value={profile}>
                        {profile}
                      </option>
                    );
                  })}
                </Control>
              </Col>
              <Col sm="3">
                <Label>Matrix</Label>
                <Control
                  as="select"
                  value={selectMatrix}
                  onChange={(e) => filterMatrix(e.target.value)}
                  custom
                >
                  <option value="0" disabled>
                    Select
                  </option>
                  {matrixOptions.map((matrix, index) => {
                    return (
                      <option key={index} value={matrix}>
                        {matrix}
                      </option>
                    );
                  })}
                </Control>
              </Col>
              <Col sm="3">
                <Label>Tag</Label>
                <Control
                  as="select"
                  value={selectTag}
                  onChange={(e) => filterTag(e.target.value)}
                  custom
                >
                  <option value="0" disabled>
                    Select
                  </option>
                  {tagOptions.map((tag, index) => {
                    return (
                      <option key={index} value={tag}>
                        {tag}
                      </option>
                    );
                  })}
                </Control>
              </Col>
            </Row>
          </div>
        </Group>
      </Form>
      <div className="d-flex">
        <a className="ml-2 px-2 py-1" href={pyPlotURL} download={getPlotName()}>
          Download Plot
        </a>
        <span className="ml-auto">
          <LoadingOverlay active={downloadOverlay} />
          <Button
            className="px-2 py-1"
            variant="link"
            onClick={() => downloadResults()}
          >
            Download Results
          </Button>
        </span>
      </div>
      <div className="border rounded p-2 mb-2">
        <Row>
          <Col>
            <img className="w-100 my-4" src={pyPlotURL}></img>
          </Col>
        </Row>
      </div>
      <div className="border rounded p-1 my-2">
        <div>python</div>
        <div>stdout</div>
        <pre className="border">{debugPy.stdout}</pre>
        <div>stderr</div>
        <pre className="border">{debugPy.stderr}</pre>
      </div>
    </div>
  );
}
