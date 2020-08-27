import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from 'react-select';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchMutationalProfiles,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;

export default function MutationalProfiles() {
  const [downloadOverlay, setOverlay] = useState(false);
  const { svgList, displayTab } = useSelector(
    (state) => state.visualizeResults
  );
  const rootURL = window.location.pathname;
  const {
    filtered,
    selectName,
    selectProfile,
    selectMatrix,
    selectFilter,
    nameOptions,
    profileOptions,
    matrixOptions,
    filterOptions,
    plotURL,
    debug,
    displayDebug,
  } = useSelector((state) => state.mutationalProfiles);

  // set inital plot
  useEffect(() => {
    if (!plotURL && displayTab == 'mutationalProfiles') {
      setPyPlot();
    }
  }, [filtered]);
  // set new plots on dropdown change
  useEffect(() => {
    if (plotURL && displayTab == 'mutationalProfiles') {
      setPyPlot();
    }
  }, [selectName, selectProfile, selectMatrix, selectFilter]);

  function getPlotName() {
    if (filtered.length) {
      const plot = filtered[0];

      return `${plot.Sample_Name}-${plot.Profile_Type}-${plot.Matrix_Size}-${plot.Filter}.svg`;
    } else return '';
  }

  async function setPyPlot() {
    if (filtered.length) {
      const plot = filtered[0];
      try {
        const response = await fetch(`${rootURL}visualize/svg`, {
          method: 'POST',
          headers: {
            Accept: 'image/svg',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ path: plot.Path }),
        });
        if (!response.ok) {
          const msg = await response.text();
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
    const filteredPlots = svgList.filter((plot) => plot.Sample_Name == name);
    const profileOptions = [
      ...new Set(filteredPlots.map((plot) => plot.Profile_Type)),
    ];
    const matrixOptions = [
      ...new Set(
        filteredPlots
          .filter((plot) => plot.Profile_Type == profileOptions[0])
          .map((plot) => plot.Matrix_Size)
      ),
    ];
    const filterOptions = [
      ...new Set(
        filteredPlots
          .filter((plot) => plot.Matrix_Size == matrixOptions[0])
          .map((plot) => plot.Filter)
      ),
    ];

    dispatchMutationalProfiles({
      selectName: name,
      selectProfile: profileOptions[0],
      selectMatrix: matrixOptions[0],
      selectFilter: filterOptions[0],
      profileOptions: profileOptions,
      matrixOptions: matrixOptions,
      filterOptions: filterOptions,
      filtered: filteredPlots,
    });
  }

  function filterProfileType(profile) {
    const filteredPlots = svgList.filter(
      (plot) => plot.Sample_Name == selectName && plot.Profile_Type == profile
    );
    const matrixOptions = [
      ...new Set(filteredPlots.map((plot) => plot.Matrix_Size)),
    ];
    const filterOptions = [
      ...new Set(
        filteredPlots
          .filter((plot) => plot.Matrix_Size == matrixOptions[0])
          .map((plot) => plot.Filter)
      ),
    ];

    dispatchMutationalProfiles({
      selectProfile: profile,
      selectMatrix: matrixOptions[0],
      selectFilter: filterOptions[0],
      matrixOptions: matrixOptions,
      filterOptions: filterOptions,
      filtered: filteredPlots,
    });
  }

  function filterMatrix(matrix) {
    const filteredPlots = svgList.filter(
      (plot) =>
        plot.Sample_Name == selectName &&
        plot.Profile_Type == selectProfile &&
        plot.Matrix_Size == matrix
    );
    const filterOptions = [
      ...new Set(filteredPlots.map((plot) => plot.Filter)),
    ];

    dispatchMutationalProfiles({
      selectMatrix: matrix,
      selectFilter: filterOptions[0],
      filterOptions: filterOptions,
      filtered: filteredPlots,
    });
  }

  function filterTag(tag) {
    const filteredPlots = svgList.filter(
      (plot) =>
        plot.Sample_Name == selectName &&
        plot.Profile_Type == selectProfile &&
        plot.Matrix_Size == selectMatrix &&
        plot.Filter == tag
    );

    dispatchMutationalProfiles({
      selectFilter: tag,
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
                <Select
                  options={nameOptions}
                  value={[selectName]}
                  onChange={(name) => {
                    dispatchMutationalProfiles({ selectName: name || '' });
                    if (nameOptions.indexOf(name) > -1) filterSampleName(name);
                  }}
                  getOptionLabel={(option) => option}
                  getOptionValue={(option) => option}
                />
              </Col>
              <Col sm="3">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={selectProfile}
                  onChange={(e) => filterProfileType(e.target.value)}
                  custom
                  disabled={!selectName}
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
                  disabled={!selectName || !selectProfile}
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
                <Label>Filter</Label>
                <Control
                  as="select"
                  value={selectFilter}
                  onChange={(e) => filterTag(e.target.value)}
                  custom
                  disabled={!selectName || !selectProfile || !selectMatrix}
                >
                  <option value="0" disabled>
                    Select
                  </option>
                  {filterOptions.map((tag, index) => {
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
            <img className="w-100 my-4 h-500" src={plotURL}></img>
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
