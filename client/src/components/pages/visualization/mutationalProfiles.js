import React, { useEffect } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import Plot from '../../controls/plot/plot';
import Select from '../../controls/select/select';
import {
  value2d,
  filter2d,
  unique2d,
  defaultProfile,
  defaultMatrix,
  defaultFilter,
} from '../../../services/utils';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';

const actions = { ...visualizationActions, ...modalActions };

export default function MutationalProfiles() {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeMutationalProfiles = (state) =>
    dispatch(actions.mergeVisualization({ mutationalProfiles: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { source } = visualization.visualize;
  const { projectID, svgList, displayTab } = visualization.results;
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
  } = visualization.mutationalProfiles;

  const { columns } = svgList;

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

  function getdownloadName() {
    if (filtered.length) {
      const plot = filtered[0];
      if (source == 'user')
        return `${value2d(plot, 'Sample_Name', columns)}-${value2d(
          'Profile_Type',
          plot,
          columns
        )}-${value2d(plot, 'Matrix_Size', columns)}-${value2d(
          'Filter',
          plot,
          columns
        )}.svg`;
      else return value2d(plot, 'Sample', columns).split('/').slice(-1)[0];
    } else return '';
  }

  async function setPlot() {
    if (filtered.length) {
      const plot = filtered[0];

      try {
        const response =
          source == 'user'
            ? await fetch(
                `api/results/${projectID}/${value2d(plot, 'Path', columns)}`
              )
            : // : await fetch(`api/public/${value2d(plot, 'Path', columns)}`);
              await fetch(`api/getImageS3`, {
                method: 'POST',
                headers: {
                  Accept: 'image/svg',
                  'Content-Type': 'application/json',
                },
                body: JSON.stringify({ path: value2d(plot, 'Path', columns) }),
              });

        if (!response.ok) {
          const msg = await response.text();
          mergeError(msg);
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL.length) URL.revokeObjectURL(plotURL);

          mergeMutationalProfiles({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
    }
  }

  function handleSample(name) {
    if (source == 'user') {
      const filteredPlots = filter2d(name, svgList.data);
      const profileOptions = unique2d('Profile_Type', columns, filteredPlots);
      const profile = defaultProfile(profileOptions);
      const matrixOptions = unique2d(
        'Matrix_Size',
        columns,
        filter2d(profile, filteredPlots)
      );
      const matrix = defaultMatrix(profile, matrixOptions);
      const filterOptions = unique2d(
        'Filter',
        columns,
        filter2d(matrix, filteredPlots)
      );

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
      const filteredPlots = filter2d(name, svgList.data);
      const profileOptions = [
        ...new Set(
          filteredPlots.map(
            (row) => value2d(row, 'Profile', columns).match(/[a-z]+/gi)[0]
          )
        ),
      ];
      const profile = defaultProfile(profileOptions);
      const matrixOptions = [
        ...new Set(
          filteredPlots
            .filter((row) => value2d(row, 'Profile', columns).includes(profile))
            .map((row) => value2d(row, 'Profile', columns).match(/\d+/gi)[0])
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
      const filteredPlots = filter2d([selectName, profile], svgList.data);
      const matrixOptions = unique2d('Matrix_Size', columns, filteredPlots);
      const matrix = defaultMatrix(profile, matrixOptions);
      const filterOptions = unique2d(
        'Filter',
        columns,
        filter2d(matrix, filteredPlots)
      );

      mergeMutationalProfiles({
        selectProfile: profile,
        selectMatrix: matrix,
        selectFilter: defaultFilter(filterOptions),
        matrixOptions: matrixOptions,
        filterOptions: filterOptions,
        filtered: filteredPlots,
      });
    } else {
      const filteredPlots = svgList.data.filter(
        (row) =>
          row.includes(selectName) &&
          value2d(row, 'Profile', columns).indexOf(profile) > -1
      );
      const matrixOptions = [
        ...new Set(
          filteredPlots.map(
            (row) => value2d(row, 'Profile', columns).match(/\d+/gi)[0]
          )
        ),
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
      const filteredPlots = filter2d(
        [selectName, selectProfile, matrix],
        svgList.data
      );
      const filterOptions = unique2d('Filter', columns, filteredPlots);

      mergeMutationalProfiles({
        selectMatrix: matrix,
        selectFilter: defaultFilter(filterOptions),
        filterOptions: filterOptions,
        filtered: filteredPlots,
      });
    } else {
      const filteredPlots = filter2d(
        [selectName, `${selectProfile + matrix}`],
        svgList.data
      );

      mergeMutationalProfiles({
        selectMatrix: matrix,
        filtered: filteredPlots,
      });
    }
  }

  function handleTag(tag) {
    const filteredPlots = filter2d(
      [selectName, selectProfile, selectMatrix, tag],
      svgList.data
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
          <Col lg="3">
            <Select
              id="mpSampleName"
              label="Sample Name"
              value={selectName}
              options={nameOptions}
              onChange={handleSample}
            />
          </Col>
          <Col lg="3">
            <Select
              id="mpProfileType"
              label="Profile Type"
              value={selectProfile}
              options={profileOptions}
              onChange={handleProfile}
            />
          </Col>
          <Col lg="3">
            <Select
              id="mpMatrixSize"
              label="Matrix Size"
              value={selectMatrix}
              options={matrixOptions}
              onChange={handleMatrix}
            />
          </Col>
          <Col lg="3">
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
      <Plot
        className="p-3"
        downloadName={getdownloadName()}
        plotURL={plotURL}
      />
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
