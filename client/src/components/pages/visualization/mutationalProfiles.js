import React, { useEffect, useState } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import Plot from '../../controls/plot/plot';
import CustomSelect from '../../controls/select/select-old';
import Description from '../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import {
  defaultProfile,
  defaultMatrix,
  defaultFilter,
} from '../../../services/utils';
import { NavHashLink } from 'react-router-hash-link';

const actions = { ...visualizationActions, ...modalActions };

export default function MutationalProfiles() {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeMutationalProfiles = (state) =>
    dispatch(actions.mergeVisualization({ mutationalProfiles: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { source, projectID, svgList, displayTab } = visualization.main;
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
    plotPath,
  } = visualization.mutationalProfiles;

  const [loading, setLoading] = useState(false);

  // populate controls
  useEffect(() => {
    if (!nameOptions.length) {
      const names = [
        ...new Set(
          svgList
            .filter((row) => row.Path != 'msigportal/Database/Seqmatrix/NA')
            .map((row) => row.Sample || row.Sample_Name)
        ),
      ];

      mergeMutationalProfiles({ nameOptions: names });
      handleSample(names[0]);
    }
  }, [selectName]);
  // set inital plot
  useEffect(() => {
    if (!plotPath && displayTab == 'mutationalProfiles') {
      setPlot();
    }
  }, [filtered, displayTab]);
  // set new plots on dropdown change
  useEffect(() => {
    if (plotPath && displayTab == 'mutationalProfiles') {
      setPlot();
    }
  }, [selectName, selectProfile, selectMatrix, selectFilter]);

  function getdownloadName() {
    if (filtered.length) {
      const plot = filtered[0];
      if (source == 'user')
        return `${plot.Sample_Name}-${plot.Profile_Type}-${plot.Matrix_Size}-${plot.Filter}.svg`;
      else return plot.Sample.split('/').slice(-1)[0];
    } else return '';
  }

  async function setPlot() {
    if (filtered.length) {
      setLoading(true);
      const plot = filtered[0];

      try {
        const response =
          source == 'user'
            ? await fetch(`web/results/${projectID}${plot.Path}`)
            : await fetch(`web/getImageS3`, {
                method: 'POST',
                headers: {
                  Accept: 'image/svg',
                  'Content-Type': 'application/json',
                },
                body: JSON.stringify({ path: plot.Path }),
              });

        if (!response.ok) {
          const msg = await response.text();
          // console.log(msg);
          // mergeError(msg);
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotPath.length) URL.revokeObjectURL(plotPath);

          mergeMutationalProfiles({
            plotPath: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
      setLoading(false);
    }
  }

  function handleSample(name) {
    if (source == 'user') {
      const filteredPlots = svgList.filter((row) => row.Sample_Name == name);
      const profileOptions = [
        ...new Set(filteredPlots.map(({ Profile_Type }) => Profile_Type)),
      ];
      const profile = defaultProfile(profileOptions);
      const matrixOptions = [
        ...new Set(
          filteredPlots
            .filter((row) => row.Profile_Type == profile)
            .map((row) => row.Matrix_Size)
        ),
      ].sort((a, b) => a - b);
      const matrix = defaultMatrix(profile, matrixOptions);
      const filterOptions = [
        ...new Set(
          filteredPlots
            .filter(
              (row) => row.Matrix_Size == matrix && row.Profile_Type == profile
            )
            .map((row) => row.Filter)
        ),
      ];

      mergeMutationalProfiles({
        selectName: name,
        selectProfile: profile,
        selectMatrix: matrix,
        selectFilter: defaultFilter(filterOptions),
        profileOptions: profileOptions,
        matrixOptions: matrixOptions,
        filterOptions: filterOptions,
        filtered: filteredPlots,
      });
    } else {
      const filteredPlots = svgList.filter(
        (row) =>
          row.Sample == name && row.Path != 'msigportal/Database/Seqmatrix/NA'
      );
      const profileOptions = [
        ...new Set(
          filteredPlots.map((row) => row.Profile.match(/[a-z]+/gi)[0])
        ),
      ];
      const profile = defaultProfile(profileOptions);
      const matrixOptions = [
        ...new Set(
          filteredPlots
            .filter((row) => row.Profile.includes(profile))
            .map((row) => row.Profile.match(/\d+/gi)[0])
        ),
      ].sort((a, b) => a - b);

      const matrix = defaultMatrix(profile, matrixOptions);

      mergeMutationalProfiles({
        selectName: name,
        selectProfile: profile,
        selectMatrix: matrix,
        profileOptions: profileOptions,
        matrixOptions: matrixOptions,
        filtered: filteredPlots,
      });
    }
  }

  function handleProfile(profile) {
    if (source == 'user') {
      const filteredPlots = svgList.filter(
        (row) => row.Sample_Name == selectName && row.Profile_Type == profile
      );
      const matrixOptions = [
        ...new Set(filteredPlots.map((row) => row.Matrix_Size)),
      ].sort((a, b) => a - b);
      const matrix = defaultMatrix(profile, matrixOptions);
      const filterOptions = [
        ...new Set(
          filteredPlots
            .filter((row) => row.Matrix_Size == matrix)
            .map((row) => row.Filter)
        ),
      ];

      mergeMutationalProfiles({
        selectProfile: profile,
        selectMatrix: matrix,
        selectFilter: defaultFilter(filterOptions),
        matrixOptions: matrixOptions,
        filterOptions: filterOptions,
        filtered: filteredPlots,
      });
    } else {
      const filteredPlots = svgList.filter(
        (row) =>
          row.Sample == selectName &&
          row.Profile.indexOf(profile) > -1 &&
          row.Path != 'msigportal/Database/Seqmatrix/NA'
      );
      const matrixOptions = [
        ...new Set(filteredPlots.map((row) => row.Profile.match(/\d+/gi)[0])),
      ].sort((a, b) => a - b);

      mergeMutationalProfiles({
        selectProfile: profile,
        selectMatrix: defaultMatrix(profile, matrixOptions),
        matrixOptions: matrixOptions,
        filtered: filteredPlots,
      });
    }
  }

  function handleMatrix(matrix) {
    if (source == 'user') {
      const filteredPlots = svgList.filter(
        (row) =>
          row.Sample_Name == selectName &&
          row.Profile_Type == selectProfile &&
          row.Matrix_Size == matrix
      );
      const filterOptions = [
        ...new Set(filteredPlots.map((row) => row.Filter)),
      ];

      mergeMutationalProfiles({
        selectMatrix: matrix,
        selectFilter: defaultFilter(filterOptions),
        filterOptions: filterOptions,
        filtered: filteredPlots,
      });
    } else {
      const filteredPlots = svgList.filter(
        (row) =>
          row.Sample == selectName &&
          row.Profile == selectProfile + matrix &&
          row.Path != 'msigportal/Database/Seqmatrix/NA'
      );

      mergeMutationalProfiles({
        selectMatrix: matrix,
        filtered: filteredPlots,
      });
    }
  }

  function handleTag(tag) {
    const filteredPlots = svgList.filter(
      (row) =>
        row.Sample_Name == selectName &&
        row.Profile_Type == selectProfile &&
        row.Matrix_Size == selectMatrix &&
        row.Filter == tag
    );

    mergeMutationalProfiles({
      selectFilter: tag,
      filtered: filteredPlots,
    });
  }

  return (
    <div className="bg-white border rounded">
      <div className="p-3">
        <Description
          className="m-0"
          less={
            <span>
              Below you can visualize different mutational profiles for a given
              sample.
            </span>
          }
          more={
            <span>
              {' '}
              Use the dropdown arrow to select the [Sample Name], [Profile Type]
              and [Matrix Size]. The [Filter] option is only available if [Split
              Mutations According to Filter] is selected while analyzing user
              data. For additional information on [Profile Type] and [Matrix
              Size], click <NavHashLink to="/faq#sbs">here</NavHashLink>.
            </span>
          }
        />
      </div>
      <hr />
      <Form className="p-3">
        <Row>
          <Col lg="auto">
            <CustomSelect
              id="mpSampleName"
              label="Sample Name"
              value={selectName}
              options={nameOptions}
              onChange={handleSample}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mpProfileType"
              label="Profile Type"
              value={selectProfile}
              options={profileOptions}
              onChange={handleProfile}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mpMatrixSize"
              label="Matrix Size"
              value={selectMatrix}
              options={matrixOptions}
              onChange={handleMatrix}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mpFilter"
              label="Filter"
              value={selectFilter}
              options={filterOptions}
              onChange={handleTag}
              disabled={source == 'public'}
            />
          </Col>
        </Row>
      </Form>

      <hr />
      <div>
        <LoadingOverlay active={loading} />
        <Plot
          className="p-3"
          downloadName={getdownloadName()}
          plotPath={plotPath}
          height="500px"
        />
      </div>
    </div>
  );
}
