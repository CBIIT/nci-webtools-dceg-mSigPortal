import React, { useEffect, useState } from 'react';
import { Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualize,
  dispatchVisualizeResults,
  dispatchMutationalProfiles,
  dispatchCosineSimilarity,
  dispatchProfileComparison,
  dispatchPCA,
} from '../../../services/store';

import ProfilerSummary from './profilerSummary';
import MutationalProfiles from './mutationalProfiles';
import CosineSimilarity from './cosineSimilarity';
import MutationalPattern from './mutationalPattern';
import ProfileComparison from './profileComparison';
import PCA from './pca';
import Download from './download';

const { Header, Body } = Card;
const { Item, Link } = Nav;

export default function Results({ setOpenSidebar }) {
  const { error, projectID, displayTab, svgList, matrixList } = useSelector(
    (state) => state.visualizeResults
  );
  const { source } = useSelector((state) => state.visualize);
  const mutationalProfiles = useSelector((state) => state.mutationalProfiles);
  const { signatureSetOptions } = useSelector((state) => state.pca);
  const rootURL = window.location.pathname;
  const [retrieveSvgList, setAttempt] = useState(false);
  // get mapping of plots after retrieving projectID
  useEffect(() => {
    if (source == 'user') {
      if (
        projectID &&
        !retrieveSvgList &&
        !mutationalProfiles.filtered.length
      ) {
        setAttempt(true);
        getSummary();
      } else if (projectID && !signatureSetOptions.length) mapSvgList();
    } else {
      if (svgList.length > 0 && !mutationalProfiles.filtered.length)
        mapPublicData();
    }
  }, [svgList, projectID, source]);

  // reload summary information
  async function getSummary() {
    console.log(projectID);
    const response = await fetch(`${rootURL}getSummary`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ projectID: projectID }),
    });
    const {
      svgList,
      statistics,
      matrixList,
      downloads,
    } = await response.json();
    dispatchVisualizeResults({
      svgList: svgList,
      statistics: statistics,
      matrixList: matrixList,
      downloads: downloads,
    });
  }

  // retrieve mapping of samples to plots from svgList file
  async function mapSvgList() {
    dispatchVisualize({
      loading: {
        active: true,
      },
    });

    const nameOptions = [...new Set(svgList.map((plot) => plot.Sample_Name))];
    const profileOptions = [
      ...new Set(svgList.map((plot) => plot.Profile_Type)),
    ];
    const matrixOptions = [...new Set(svgList.map((plot) => plot.Matrix_Size))];
    const filterOptions = [...new Set(svgList.map((plot) => plot.Filter))];

    const selectName = mutationalProfiles.selectName || nameOptions[0];
    const selectProfile = mutationalProfiles.selectProfile || profileOptions[0];
    const selectMatrix = mutationalProfiles.selectMatrix || matrixOptions[0];
    const selectFilter = mutationalProfiles.selectFilter || filterOptions[0];

    const filteredPlots = svgList.filter(
      (plot) =>
        plot.Sample_Name == selectName &&
        plot.Profile_Type == selectProfile &&
        plot.Matrix_Size == selectMatrix &&
        plot.Filter == selectFilter
    );

    const filteredProfileOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Sample_Name == selectName)
          .map((plot) => plot.Profile_Type)
      ),
    ];
    const filteredMatrixOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Profile_Type == selectProfile)
          .map((plot) => plot.Matrix_Size)
      ),
    ];
    const filteredFilterOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Matrix_Size == selectMatrix)
          .map((plot) => plot.Filter)
      ),
    ];
    const filteredMatrixList = [
      ...new Set(
        matrixList
          .filter((matrix) => matrix.Profile_Type == selectProfile)
          .map((matrix) => matrix.Matrix_Size)
      ),
    ];

    const refSignatureSetOptions = await (
      await getRefSigOptions(profileOptions[0])
    ).json();

    dispatchMutationalProfiles({
      filtered: filteredPlots,
      nameOptions: nameOptions,
      profileOptions: filteredProfileOptions,
      matrixOptions: filteredMatrixOptions,
      filterOptions: filteredFilterOptions,
      selectName: selectName,
      selectProfile: selectProfile,
      selectMatrix: selectMatrix,
      selectFilter: selectFilter,
    });

    dispatchCosineSimilarity({
      withinProfileType: profileOptions[0],
      refProfileType: profileOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      withinMatrixSize: filteredMatrixList[0],
      withinMatrixOptions: filteredMatrixList,
      userProfileType: profileOptions[0],
      userMatrixSize: filteredMatrixOptions[0],
      userMatrixOptions: filteredMatrixOptions,
    });

    dispatchProfileComparison({
      withinProfileType: profileOptions[0],
      withinSampleName1: nameOptions[0],
      withinSampleName2: nameOptions[1],
      refProfileType: profileOptions[0],
      refSampleName: nameOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      userProfileType: profileOptions[0],
      userMatrixSize: filteredMatrixOptions[0],
      userMatrixOptions: filteredMatrixOptions,
      userSampleName: nameOptions[0],
    });

    dispatchPCA({
      profileType: profileOptions[0],
      signatureSet: refSignatureSetOptions[0],
      signatureSetOptions: refSignatureSetOptions,
      userProfileType: profileOptions[0],
      userMatrixSize: filteredMatrixList[0],
      userMatrixOptions: filteredMatrixOptions,
    });

    dispatchVisualize({
      loading: {
        active: false,
      },
    });
    setOpenSidebar(false);
  }

  // retrieve mapping of samples to plots from svgList file
  async function mapPublicData() {
    dispatchVisualize({
      loading: {
        active: true,
      },
    });

    const nameOptions = [...new Set(svgList.map((plot) => plot.Sample))];
    const profileOptions = [
      ...new Set(svgList.map((plot) => plot.Profile.match(/[a-z]+/gi)[0])),
    ];
    const matrixOptions = [
      ...new Set(svgList.map((plot) => plot.Profile.match(/\d+/gi)[0])),
    ];
    const filterOptions = ['NA'];

    const selectName = mutationalProfiles.selectName || nameOptions[0];
    const selectProfile = mutationalProfiles.selectProfile || profileOptions[0];
    const selectMatrix = mutationalProfiles.selectMatrix || matrixOptions[0];
    const selectFilter = filterOptions[0];

    const filteredPlots = svgList.filter(
      (plot) =>
        plot.Sample == selectName &&
        plot.Profile.indexOf(selectProfile) > -1 &&
        plot.Profile.indexOf(selectMatrix) > -1
    );

    const filteredProfileOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Sample == selectName)
          .map((plot) => plot.Profile.match(/[a-z]+/gi)[0])
      ),
    ];

    const filteredMatrixOptions = [
      ...new Set(
        svgList
          .filter(
            (plot) =>
              plot.Sample == selectName &&
              plot.Profile.indexOf(filteredProfileOptions[0]) > -1
          )
          .map((plot) => plot.Profile.match(/\d+/gi)[0])
      ),
    ];

    const refSignatureSetOptions = await (
      await getRefSigOptions(profileOptions[0])
    ).json();

    dispatchMutationalProfiles({
      filtered: filteredPlots,
      nameOptions: nameOptions,
      profileOptions: filteredProfileOptions,
      matrixOptions: filteredMatrixOptions,
      filterOptions: filterOptions,
      selectName: selectName,
      selectProfile: selectProfile,
      selectMatrix: selectMatrix,
      selectFilter: selectFilter,
    });

    dispatchCosineSimilarity({
      withinProfileType: profileOptions[0],
      refProfileType: profileOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      withinMatrixSize: filteredMatrixOptions[0],
      withinMatrixOptions: filteredMatrixOptions,
    });

    dispatchProfileComparison({
      withinProfileType: profileOptions[0],
      withinSampleName1: nameOptions[0],
      withinSampleName2: nameOptions[1],
      refProfileType: profileOptions[0],
      refSampleName: nameOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
    });

    dispatchPCA({
      profileType: profileOptions[0],
      signatureSet: refSignatureSetOptions[0],
      signatureSetOptions: refSignatureSetOptions,
    });

    dispatchVisualize({
      loading: {
        active: false,
      },
    });
    setOpenSidebar(false);
  }

  function submitR(fn, args) {
    return fetch(`${rootURL}visualizeR`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ fn: fn, args: args, projectID: projectID }),
    });
  }

  function getRefSigOptions(profileType) {
    return fetch(`${rootURL}visualizeR/getReferenceSignatureSets`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ profileType: profileType }),
    });
  }

  //   download text results files
  async function downloadResults(txtPath) {
    try {
      const response = await fetch(`${rootURL}downloadPlotData`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ path: txtPath }),
      });
      if (!response.ok) {
        const { msg } = await response.json();
        dispatchError(msg);
      } else {
        const file = await response.blob();
        const objectURL = URL.createObjectURL(file);
        const tempLink = document.createElement('a');

        tempLink.href = objectURL;
        tempLink.download = txtPath.split('/').slice(-1)[0];
        document.body.appendChild(tempLink);
        tempLink.click();
        document.body.removeChild(tempLink);
        URL.revokeObjectURL(objectURL);
      }
    } catch (err) {
      dispatchError(err);
    }
  }

  return error.length ? (
    <h4 className="text-danger">{error}</h4>
  ) : mutationalProfiles.filtered.length ? (
    <Card>
      <Header>
        <Nav variant="pills" defaultActiveKey="#mutationalProfiles">
          {[
            { title: 'Profiler Summary', id: 'profilerSummary' },
            { title: 'Mutational Profiles', id: 'mutationalProfiles' },
            { title: 'Cosine Similarity', id: 'cosineSimilarity' },
            {
              title: 'Mutational Pattern Enrichment Analysis',
              id: 'mutationalPattern',
            },
            { title: 'Profile Comparison', id: 'profileComparison' },
            { title: 'PCA', id: 'pca' },
            { title: 'Download', id: 'download' },
          ].map(({ title, id }) => {
            return (
              <Item key={id}>
                <Link
                  active={displayTab == id}
                  onClick={() => dispatchVisualizeResults({ displayTab: id })}
                >
                  {title}
                </Link>
              </Item>
            );
          })}
        </Nav>
      </Header>
      <Body
        style={{
          display: displayTab == 'profilerSummary' ? 'block' : 'none',
        }}
      >
        <ProfilerSummary submitR={(fn, args) => submitR(fn, args)} />
      </Body>
      <Body
        style={{
          display: displayTab == 'mutationalProfiles' ? 'block' : 'none',
        }}
      >
        <MutationalProfiles />
      </Body>
      <Body
        style={{ display: displayTab == 'cosineSimilarity' ? 'block' : 'none' }}
      >
        <CosineSimilarity
          getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
          downloadResults={(path) => downloadResults(path)}
          submitR={(fn, args) => submitR(fn, args)}
        />
      </Body>
      <Body
        style={{
          display: displayTab == 'mutationalPattern' ? 'block' : 'none',
        }}
      >
        <MutationalPattern
          downloadResults={(path) => downloadResults(path)}
          submitR={(fn, args) => submitR(fn, args)}
        />
      </Body>

      <Body
        style={{
          display: displayTab == 'profileComparison' ? 'block' : 'none',
        }}
      >
        <ProfileComparison
          getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
          submitR={(fn, args) => submitR(fn, args)}
        />
      </Body>
      <Body style={{ display: displayTab == 'pca' ? 'block' : 'none' }}>
        <PCA
          getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
          downloadResults={(path) => downloadResults(path)}
          submitR={(fn, args) => submitR(fn, args)}
        />
      </Body>
      <Body style={{ display: displayTab == 'download' ? 'block' : 'none' }}>
        <Download />
      </Body>
    </Card>
  ) : (
    <div>
      <h2>Instructions</h2>
      <p>Upload Sample Variants and specify parameters</p>
    </div>
  );
}
