import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchMutationalProfiles,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;

export default function MutationalProfiles() {
  const [downloadOverlay, setOverlay] = useState(false);
  const { summary } = useSelector((state) => state.visualizeResults);
  const rootURL = window.location.pathname;
  const {
    filtered,
    selectName,
    selectProfile,
    selectMatrix,
    selectTag,
    nameOptions,
    profileOptions,
    matrixOptions,
    tagOptions,
    plotURL,
    debug,
    displayDebug,
  } = useSelector((state) => state.mutationalProfiles);

  // set inital plot
  useEffect(() => {
    if (!plotURL.length) {
      setPyPlot();
    }
  }, [filtered]);

  // change plot
  useEffect(() => {
    if (filtered.length && plotURL.length) {
      setPyPlot();
    }
  }, [selectName, selectProfile, selectMatrix, selectTag]);

  function getPlotName() {
    if (filtered.length) {
      const plot = filtered[0];

      return `${plot.Sample_Name}-${plot.Profile_Type}-${plot.Matrix}-${plot.Tag}.svg`;
    } else return '';
  }

  async function setPyPlot() {
    if (filtered.length) {
      const plot = filtered[0];
      try {
        const response = await fetch(`${rootURL}visualize/svg`, {
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

          dispatchMutationalProfiles({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      console.log('filtered summary empty');
    }
  }

  // async function downloadResults() {
  //   setOverlay(true);
  //   try {
  //     const response = await fetch(`${rootURL}visualize/results`, {
  //       method: 'POST',
  //       headers: {
  //         Accept: 'application/json',
  //         'Content-Type': 'application/json',
  //       },
  //       body: JSON.stringify({ projectID: projectID }),
  //     });
  //     if (!response.ok) {
  //       setOverlay(false);
  //       const { msg } = await response.json();
  //       dispatchError(msg);
  //     } else {
  //       setOverlay(false);
  //       const req = await response.json();
  //       const id = req.projectID;
  //       const tempLink = document.createElement('a');

  //       tempLink.href = `/visualize/download?id=${id}`;
  //       document.body.appendChild(tempLink);
  //       tempLink.click();
  //       document.body.removeChild(tempLink);
  //     }
  //   } catch (err) {
  //     dispatchError(err);
  //   }
  // }

  function filterSampleName(name) {
    const filteredPlots = summary.filter((plot) => plot.Sample_Name == name);
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

    dispatchMutationalProfiles({
      selectName: name,
      selectProfile: profileOptions[0],
      selectMatrix: matrixOptions[0],
      selectTag: tagOptions[0],
      profileOptions: profileOptions,
      matrixOptions: matrixOptions,
      tagOptions: tagOptions,
      filtered: filteredPlots,
    });
  }

  function filterProfileType(profile) {
    const filteredPlots = summary.filter(
      (plot) => plot.Sample_Name == selectName && plot.Profile_Type == profile
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

    dispatchMutationalProfiles({
      selectProfile: profile,
      selectMatrix: matrixOptions[0],
      selectTag: tagOptions[0],
      matrixOptions: matrixOptions,
      tagOptions: tagOptions,
      filtered: filteredPlots,
    });
  }

  function filterMatrix(matrix) {
    const filteredPlots = summary.filter(
      (plot) =>
        plot.Sample_Name == selectName &&
        plot.Profile_Type == selectProfile &&
        plot.Matrix == matrix
    );
    const tagOptions = [...new Set(filteredPlots.map((plot) => plot.Tag))];

    dispatchMutationalProfiles({
      selectMatrix: matrix,
      selectTag: tagOptions[0],
      tagOptions: tagOptions,
      filtered: filteredPlots,
    });
  }

  function filterTag(tag) {
    const filteredPlots = summary.filter(
      (plot) =>
        plot.Sample_Name == selectName &&
        plot.Profile_Type == selectProfile &&
        plot.Matrix == selectMatrix &&
        plot.Tag == tag
    );

    dispatchMutationalProfiles({
      selectTag: tag,
      filtered: filteredPlots,
    });
  }

  return (
    <div>
      <Form>
        <Group controlId="selectPlot">
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
                <Label>Matrix Size</Label>
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
        <a className="px-2 py-1" href={plotURL} download={getPlotName()}>
          Download Plot
        </a>
        <span className="ml-auto">
          <LoadingOverlay active={downloadOverlay} />
        </span>
      </div>
      <div className="border rounded p-2 mb-2">
        <Row>
          <Col>
            <img className="w-100 my-4" src={plotURL}></img>
          </Col>
        </Row>
      </div>
      <Button
        variant="link"
        className="p-0 mt-5"
        onClick={() =>
          dispatchMutationalProfiles({
            displayDebug: !displayDebug,
          })
        }
      >
        Python Debug
      </Button>
      <pre
        className="border rounded p-1 "
        style={{ display: displayDebug ? 'block' : 'none' }}
      >
        <div>stdout</div>
        <pre className="border">{debug.stdout}</pre>
        <div>stderr</div>
        <pre className="border">{debug.stderr}</pre>
      </pre>
    </div>
  );
}
