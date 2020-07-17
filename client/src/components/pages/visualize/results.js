import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button, Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualizeResults,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;
const { Header, Body } = Card;
const { Item, Link } = Nav;

const root =
  process.env.NODE_ENV === 'development'
    ? 'http://localhost:8330/'
    : window.location.pathname;

export default function Results() {
  const [downloadOverlay, setOverlay] = useState(false);
  const {
    projectID,
    displayTab,
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
    displayedPlotIndex,
    plotURL,
    error,
    sigProfileType,
    signatureSet,
    selectName2,
    selectSigFormula,
    sigFormula,
    rPlots,
    rPlotIndex,
    rPlotURL,
    submitOverlay,
    debug,
    debugR,
  } = useSelector((state) => state.visualizeResults);

  // get mapping after retrieving projectID
  useEffect(() => {
    if (projectID.length) {
      getSummary(projectID);
    }
  }, [projectID]);

  // load first plots after mapping is loaded
  useEffect(() => {
    if (filtered.length) {
      setPlot(0);
      submitR();
    }
  }, [mapping]);

  // load first r plot after they are recieved
  useEffect(() => {
    if (rPlots.length && !rPlotIndex.length) {
      setPlot(0, 'r');
    }
  }, [rPlots]);

  // set first plot after dropdown filters are changed
  useEffect(() => {
    setPlot(0);
  }, [filtered]);

  // retrieve mapping of samples to plots from summary file
  async function getSummary() {
    const response = await fetch(`${root}visualize/summary`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ projectID: projectID }),
    });
    const mapping = await response.json();
    const nameOptions = [...new Set(mapping.map((plot) => plot.Sample_Name))];
    const profileOptions = [
      ...new Set(mapping.map((plot) => plot.Profile_Type)),
    ];
    const matrixOptions = [...new Set(mapping.map((plot) => plot.Matrix))];
    const tagOptions = [...new Set(mapping.map((plot) => plot.Tag))];

    dispatchVisualizeResults({
      mapping: mapping,
      filtered: mapping,
      nameOptions: nameOptions,
      profileOptions: profileOptions,
      matrixOptions: matrixOptions,
      tagOptions: tagOptions,
      selectName: nameOptions[0],
      selectProfile: profileOptions[0],
      selectMatrix: matrixOptions[0],
      selectTag: tagOptions[0],
      selectName2: nameOptions[0],
    });
  }

  async function setPlot(index, type = 'python') {
    let plot = filtered[index];
    if (type == 'r') plot = rPlots[index];
    if (plot) {
      try {
        const response = await fetch(`${root}visualize/svg`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ path: plot.Location || plot }),
        });
        if (!response.ok) {
          const { msg } = await response.json();
          dispatchError(msg);
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL.length) URL.revokeObjectURL(plotURL);
          if (type == 'python') {
            dispatchVisualizeResults({
              displayedPlotIndex: index,
              plotURL: objectURL,
            });
          } else if (type == 'r') {
            dispatchVisualizeResults({
              rPlotIndex: index,
              rPlotURL: objectURL,
            });
          }
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      console.log('invalid index', index);
    }
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

    dispatchVisualizeResults({
      selectName: name,
      selectProfile: profileOptions[0],
      selectMatrix: '0',
      selectTag: '0',
      profileOptions: profileOptions,
      matrixOptions: [...new Set(filteredPlots.map((plot) => plot.Matrix))],
      tagOptions: [...new Set(filteredPlots.map((plot) => plot.Tag))],
      filtered: filteredPlots,
    });
  }

  function filterProfileType(profile) {
    const filteredPlots = mapping.filter(
      (plot) => plot.Sample_Name == selectName && plot.Profile_Type == profile
    );
    const matrixOptions = [
      ...new Set(filteredPlots.map((plot) => plot.Matrix)),
    ];

    dispatchVisualizeResults({
      selectProfile: profile,
      selectMatrix: matrixOptions[0],
      selectTag: '0',
      matrixOptions: matrixOptions,
      tagOptions: [...new Set(filteredPlots.map((plot) => plot.Tag))],
      filtered: filteredPlots,
    });
  }

  function filterMatrix(matrix) {
    const filteredPlots = mapping.filter(
      (plot) =>
        plot.Sample_Name == selectName &&
        plot.Profile_Type == selectProfile &&
        plot.Matrix == matrix
    );
    const tagOptions = [...new Set(filteredPlots.map((plot) => plot.Tag))];
    dispatchVisualizeResults({
      selectMatrix: matrix,
      selectTag: tagOptions[0],
      tagOptions: tagOptions,
      filtered: filteredPlots,
    });
  }

  function filterTag(tag) {
    const filteredPlots = mapping.filter(
      (plot) =>
        plot.Sample_Name == selectName &&
        plot.Profile_Type == selectProfile &&
        plot.Matrix == selectMatrix &&
        plot.Tag == tag
    );
    dispatchVisualizeResults({
      selectTag: tag,
      filtered: filteredPlots,
    });
  }

  async function submitR() {
    dispatchVisualizeResults({ submitOverlay: true });
    let args = {
      profileName: sigProfileType,
      signatureSetName: signatureSet,
      sampleName: selectName2,
      signatureName: null,
      formula: null,
      projectID: projectID,
    };
    selectSigFormula == 'signature'
      ? (args.signatureName = sigFormula)
      : (args.formula = sigFormula);

    try {
      const response = await fetch(`${root}api/visualizeR`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(args),
      });
      if (!response.ok) {
        const err = await response.text();
        console.log(err);
        dispatchVisualizeResults({
          debugR: err,
          rPlots: [],
        });
      } else {
        const data = await response.json();
        console.log(data);
        dispatchVisualizeResults({
          debugR: data.debugR,
          rPlots: data.plots,
        });
      }
    } catch (err) {
      dispatchError(err);
    } finally {
      dispatchVisualizeResults({ submitOverlay: false });
    }
  }

  return error.length ? (
    <h4 className="text-danger">{error}</h4>
  ) : mapping.length ? (
    <Card>
      <Header>
        <Nav variant="pills" defaultActiveKey="#python">
          <Item>
            <Link
              active={displayTab == 'python'}
              onClick={() => dispatchVisualizeResults({ displayTab: 'python' })}
            >
              Python
            </Link>
          </Item>
          <Item>
            <Link
              active={displayTab == 'r'}
              onClick={() => dispatchVisualizeResults({ displayTab: 'r' })}
            >
              R
            </Link>
          </Item>
        </Nav>
      </Header>
      <Body style={{ display: displayTab == 'python' ? 'block' : 'none' }}>
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
                    // defaultValue="unselected"
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
                    disabled={selectName == '0'}
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
                    disabled={selectProfile == '0'}
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
                    disabled={selectMatrix == '0'}
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
          <a className="ml-2 px-2 py-1" href={plotURL} download={getPlotName()}>
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
              <img className="w-100 my-4" src={plotURL}></img>
            </Col>
          </Row>
        </div>
        <div className="border rounded p-1 my-2">
          <div>python</div>
          <div>stdout</div>
          <pre className="border">{debug.stdout}</pre>
          <div>stderr</div>
          <pre className="border">{debug.stderr}</pre>
        </div>
      </Body>
      <Body style={{ display: displayTab == 'r' ? 'block' : 'none' }}>
        <Form>
          <Label>Additional Plots</Label>
          <div className="border rounded p-2">
            <LoadingOverlay active={submitOverlay} />
            <Row className="justify-content-center">
              <Col sm="4">
                <Group controlId="sigProfileType">
                  <Label>Signature Profile Type</Label>
                  <Control
                    as="select"
                    value={sigProfileType}
                    onChange={(e) =>
                      dispatchVisualizeResults({
                        sigProfileType: e.target.value,
                      })
                    }
                    custom
                  >
                    <option value="SBS">SBS</option>
                    <option value="DBS">DBS</option>
                    <option value="ID">ID</option>
                  </Control>
                </Group>
              </Col>
              <Col sm="4">
                <Group controlId="signatureSet">
                  <Label>Reference Set</Label>
                  <Control
                    as="select"
                    value={signatureSet}
                    onChange={(e) =>
                      dispatchVisualizeResults({ signatureSet: e.target.value })
                    }
                    custom
                  >
                    <option value="COSMIC v3 Signatures (SBS)">
                      COSMIC v3 Signatures (SBS)
                    </option>
                  </Control>
                </Group>
              </Col>
              <Col sm="4">
                <Group controlId="selectName2">
                  <Label>Sample Name</Label>
                  <Control
                    as="select"
                    value={selectName2}
                    onChange={(e) =>
                      dispatchVisualizeResults({ selectName2: e.target.value })
                    }
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
                </Group>
              </Col>
            </Row>
            <Row className="justify-content-center">
              <Col sm="6">
                <Group
                  controlId="selectSigFormula"
                  className="d-flex align-items-center"
                >
                  <Label className="mr-auto">Signature/Formula</Label>
                  <Control
                    as="select"
                    value={selectSigFormula}
                    onChange={(e) =>
                      dispatchVisualizeResults({
                        selectSigFormula: e.target.value,
                      })
                    }
                    custom
                    style={{ width: '150px' }}
                  >
                    <option value="signature">Signature</option>
                    <option value="formula">Formula</option>
                  </Control>
                </Group>
              </Col>
              <Col sm="6">
                <Control
                  type="text"
                  placeholder={selectSigFormula}
                  value={sigFormula}
                  onChange={(e) =>
                    dispatchVisualizeResults({ sigFormula: e.target.value })
                  }
                ></Control>
              </Col>
            </Row>
            <Row>
              <Col>
                <Button variant="primary" onClick={() => submitR()}>
                  Calculate
                </Button>
              </Col>
            </Row>
          </div>
        </Form>
        <div
          className="mt-2 p-2 border rounded"
          style={{ display: rPlots.length ? 'block' : 'none' }}
        >
          <Col sm="auto" className="p-0">
            <Control
              as="select"
              value={rPlotIndex}
              onChange={(e) => setPlot(e.target.value, 'r')}
              custom
            >
              {rPlots.map((plot, index) => {
                return (
                  <option key={index} value={index}>
                    {plot}
                  </option>
                );
              })}
            </Control>
          </Col>
          <Row>
            <Col>
              <img className="w-100 my-4" src={rPlotURL}></img>
            </Col>
          </Row>
        </div>
        <div className="border rounded p-1 my-2">
          <div>R output</div>
          <div className="border">
            {Array.isArray(debugR) ? (
              debugR.map((line, index) => {
                return (
                  <p className="m-0">
                    {index}| {line}
                  </p>
                );
              })
            ) : (
              <p>{debugR}</p>
            )}
          </div>
        </div>
      </Body>
    </Card>
  ) : (
    <div>
      <h2>Instructions</h2>
      <p>Upload Sample Variants and specify parameters</p>
    </div>
  );
}
