import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchMutationalProfiles,
} from '../../../services/store';
import Plot from '../../controls/plot/plot';
import Select from '../../controls/select/select';

const { Group, Label } = Form;

export default function MutationalProfiles() {
  const { source } = useSelector((state) => state.visualize);
  const { svgList, displayTab } = useSelector(
    (state) => state.visualizeResults
  );

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
      setPlot();
    }
  }, [filtered, displayTab]);
  // set new plots on dropdown change
  useEffect(() => {
    if (plotURL && displayTab == 'mutationalProfiles') {
      setPlot();
    }
  }, [selectName, selectProfile, selectMatrix, selectFilter]);

  function getPlotName() {
    if (filtered.length) {
      const plot = filtered[0];
      if (source == 'user')
        return `${plot.Sample_Name}-${plot.Profile_Type}-${plot.Matrix_Size}-${plot.Filter}.svg`;
      else return plot.Path.split('/').slice(-1)[0];
    } else return '';
  }

  async function setPlot() {
    if (filtered.length) {
      const plot = filtered[0];
      try {
        const response =
          source == 'user'
            ? await fetch(`api/results/${plot.Path}`)
            : await fetch(`api/public/${plot.Path}`);

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

  function filterSampleName(name) {
    if (source == 'user') {
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
    } else {
      const filteredPlots = svgList.filter((plot) => plot.Sample == name);
      const profileOptions = [
        ...new Set(
          filteredPlots.map((plot) => plot.Profile.match(/[a-z]+/gi)[0])
        ),
      ];
      const matrixOptions = [
        ...new Set(
          filteredPlots
            .filter((plot) => plot.Profile.indexOf(profileOptions[0]) > -1)
            .map((plot) => plot.Profile.match(/\d+/gi)[0])
        ),
      ];

      dispatchMutationalProfiles({
        selectName: name,
        selectProfile: profileOptions[0],
        selectMatrix: matrixOptions[0],
        selectFilter: filterOptions[0],
        profileOptions: profileOptions,
        matrixOptions: matrixOptions,
        filtered: filteredPlots,
      });
    }
  }

  function filterProfileType(profile) {
    if (source == 'user') {
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
    } else {
      const filteredPlots = svgList.filter(
        (plot) =>
          plot.Sample == selectName && plot.Profile.indexOf(profile) > -1
      );
      const matrixOptions = [
        ...new Set(filteredPlots.map((plot) => plot.Profile.match(/\d+/gi)[0])),
      ];
      dispatchMutationalProfiles({
        selectProfile: profile,
        selectMatrix: matrixOptions[0],
        matrixOptions: matrixOptions,
        filtered: filteredPlots,
      });
    }
  }

  function filterMatrix(matrix) {
    if (source == 'user') {
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
    } else {
      const filteredPlots = svgList.filter(
        (plot) =>
          plot.Sample == selectName &&
          plot.Profile.indexOf(selectProfile) > -1 &&
          plot.Profile.indexOf(matrix) > -1
      );

      dispatchMutationalProfiles({
        selectMatrix: matrix,
        filtered: filteredPlots,
      });
    }
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
        <div className="border rounded p-2">
          <Row className="justify-content-center">
            <Col sm="3">
              <Select
                className="mb-0"
                id="mpSampleName"
                label="Sample Name"
                value={selectName}
                options={nameOptions}
                onChange={filterSampleName}
              />
            </Col>
            <Col sm="3">
              <Select
                className="mb-0"
                id="mpProfileType"
                label="Profile Type"
                value={selectProfile}
                options={profileOptions}
                onChange={filterProfileType}
              />
            </Col>
            <Col sm="3">
              <Select
                className="mb-0"
                id="mpMatrixSize"
                label="Matrix Size"
                value={selectMatrix}
                options={matrixOptions}
                onChange={filterMatrix}
              />
            </Col>
            <Col sm="3">
              <Select
                className="mb-0"
                id="mpFilter"
                label="Filter"
                value={selectFilter}
                options={filterOptions}
                onChange={filterTag}
                disabled={source == 'public'}
              />
            </Col>
          </Row>
        </div>
      </Form>

      <Plot plotName={getPlotName()} plotURL={plotURL} />
      <Button
        variant="link"
        className="p-0 mt-5"
        onClick={() =>
          dispatchMutationalProfiles({
            displayDebug: !displayDebug,
          })
        }
      >
        Debug
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
