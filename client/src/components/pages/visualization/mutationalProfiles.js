import React, { useEffect, useState } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import Plot from '../../controls/plot/plot';
import Select from '../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import {
  defaultProfile,
  defaultMatrix,
  defaultFilter,
} from '../../../services/utils';

const actions = { ...visualizationActions, ...modalActions };

export default function MutationalProfiles() {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeMutationalProfiles = (state) =>
    dispatch(actions.mergeVisualization({ mutationalProfiles: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { source, projectID, svgList, displayTab } = visualization.state;
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
    debug,
    displayDebug,
  } = visualization.mutationalProfiles;

  const { columns } = svgList;
  const [loading, setLoading] = useState(false);

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
            ? await fetch(`api/results/${projectID}${plot.Path}`)
            : await fetch(`api/getImageS3`, {
                method: 'POST',
                headers: {
                  Accept: 'image/svg',
                  'Content-Type': 'application/json',
                },
                body: JSON.stringify({ path: plot.Path }),
              });

        if (!response.ok) {
          const msg = await response.text();
          console.log(msg);
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
      ];
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
      const filteredPlots = svgList.filter((row) => row.Sample == name);
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
      ];

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
      ];
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
        (row) => row.Sample == selectName && row.Profile.indexOf(profile) > -1
      );
      const matrixOptions = [
        ...new Set(filteredPlots.map((row) => row.Profile.match(/\d+/gi)[0])),
      ];

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
          row.Sample == selectName && row.Profile == selectProfile + matrix
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
        <p>
          Below you can generate different mutational profiles for a given
          sample. Use the dropdown arrow to select [Sample Name], [Profile Type]
          and [Matrix Size]. The [Filter] will only available when you select
          [Split Mutations According to Filter] for user data source.
        </p>
        <p>
          For additional information on [Profile Type] and [Matrix Size], click{' '}
          <a href="#faq">here.</a>
        </p>
      </div>
      <hr />
      <Form className="p-3">
        <Row>
          <Col lg="auto">
            <Select
              id="mpSampleName"
              label="Sample Name"
              value={selectName}
              options={nameOptions}
              onChange={handleSample}
            />
          </Col>
          <Col lg="auto">
            <Select
              id="mpProfileType"
              label="Profile Type"
              value={selectProfile}
              options={profileOptions}
              onChange={handleProfile}
            />
          </Col>
          <Col lg="auto">
            <Select
              id="mpMatrixSize"
              label="Matrix Size"
              value={selectMatrix}
              options={matrixOptions}
              onChange={handleMatrix}
            />
          </Col>
          <Col lg="auto">
            <Select
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
      <div style={{ minHeight: '400px' }}>
        <LoadingOverlay active={loading} />
        {plotPath && (
          <Plot
            className="p-3"
            downloadName={getdownloadName()}
            plotPath={plotPath}
            height="500px"
          />
        )}
      </div>
      {/* <Button
        variant="link"
        className="p-0 mt-5"
        onClick={() =>
          mergeMutationalProfiles({
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
      </pre> */}
    </div>
  );
}
