import React, { useEffect } from 'react';
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

import MutationalProfiles from './mutationalProfiles';
import CosineSimilarity from './cosineSimilarity';
import ProfileComparison from './profileComparison';
import PCA from './pca';
import Download from './download';

const { Header, Body } = Card;
const { Item, Link } = Nav;

export default function Results({ setOpenSidebar }) {
  const { error, projectID, displayTab, summary, matrixList } = useSelector(
    (state) => state.visualizeResults
  );
  const mutationalProfiles = useSelector((state) => state.mutationalProfiles);
  const { signatureSetOptions } = useSelector((state) => state.pca);
  const rootURL = window.location.pathname;

  // get mapping of plots after retrieving projectID
  useEffect(() => {
    if (summary.length) {
      // only set summary if signature set was not set
      if (!signatureSetOptions.length) mapSummary();
    } else {
      if (projectID.length) getSummary();
    }
  }, [summary]);

  // reload summary information
  async function getSummary() {
    console.log(projectID);
    const response = await fetch(`${rootURL}visualize/summary`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ projectID: projectID }),
    });
    const {
      summary,
      statistics,
      matrixList,
      downloads,
    } = await response.json();
    dispatchVisualizeResults({
      summary: summary,
      statistics: statistics,
      matrixList: matrixList,
      downloads: downloads,
    });
  }

  // retrieve mapping of samples to plots from summary file
  async function mapSummary() {
    dispatchVisualize({
      loading: {
        active: true,
      },
    });

    const nameOptions = [...new Set(summary.map((plot) => plot.Sample_Name))];
    const profileOptions = [
      ...new Set(summary.map((plot) => plot.Profile_Type)),
    ];
    const matrixOptions = [...new Set(summary.map((plot) => plot.Matrix))];
    const tagOptions = [...new Set(summary.map((plot) => plot.Tag))];

    const selectName = mutationalProfiles.selectName || nameOptions[0];
    const selectProfile = mutationalProfiles.selectProfile || profileOptions[0];
    const selectMatrix = mutationalProfiles.selectMatrix || matrixOptions[0];
    const selectTag = mutationalProfiles.selectTag || tagOptions[0];

    const filteredPlots = summary.filter(
      (plot) =>
        plot.Sample_Name == selectName &&
        plot.Profile_Type == selectProfile &&
        plot.Matrix == selectMatrix &&
        plot.Tag == selectTag
    );

    const filteredProfileOptions = [
      ...new Set(
        summary
          .filter((plot) => plot.Sample_Name == selectName)
          .map((plot) => plot.Profile_Type)
      ),
    ];
    const filteredMatrixOptions = [
      ...new Set(
        summary
          .filter((plot) => plot.Profile_Type == selectProfile)
          .map((plot) => plot.Matrix)
      ),
    ];
    const filteredTagOptions = [
      ...new Set(
        summary
          .filter((plot) => plot.Matrix == selectMatrix)
          .map((plot) => plot.Tag)
      ),
    ];
    const filteredMatrixList = [
      ...new Set(
        matrixList
          .filter((matrix) => matrix.Catalog == selectProfile)
          .map((matrix) => matrix.Type)
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
      tagOptions: filteredTagOptions,
      selectName: selectName,
      selectProfile: selectProfile,
      selectMatrix: selectMatrix,
      selectTag: selectTag,
    });

    dispatchCosineSimilarity({
      withinProfileType: profileOptions[0],
      refProfileType: profileOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      withinMatrixSize: filteredMatrixList[0],
      withinMatrixOptions: filteredMatrixList,
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
    return fetch(`${rootURL}api/visualizeR`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ fn: fn, args: args, projectID: projectID }),
    });
  }

  function getRefSigOptions(profileType) {
    return fetch(`${rootURL}api/visualizeR/getReferenceSignatureSets`, {
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
      const response = await fetch(`${rootURL}visualize/txt`, {
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
          <Item>
            <Link
              active={displayTab == 'mutationalProfiles'}
              onClick={() =>
                dispatchVisualizeResults({ displayTab: 'mutationalProfiles' })
              }
            >
              Mutational Profiles
            </Link>
          </Item>
          <Item>
            <Link
              active={displayTab == 'cosineSimilarity'}
              onClick={() =>
                dispatchVisualizeResults({ displayTab: 'cosineSimilarity' })
              }
            >
              Cosine Similarity
            </Link>
          </Item>
          <Item>
            <Link
              active={displayTab == 'profileComparison'}
              onClick={() =>
                dispatchVisualizeResults({ displayTab: 'profileComparison' })
              }
            >
              Profile Comparison
            </Link>
          </Item>
          <Item>
            <Link
              active={displayTab == 'pca'}
              onClick={() => dispatchVisualizeResults({ displayTab: 'pca' })}
            >
              PCA
            </Link>
          </Item>
          <Item>
            <Link
              active={displayTab == 'download'}
              onClick={() =>
                dispatchVisualizeResults({ displayTab: 'download' })
              }
            >
              Download
            </Link>
          </Item>
        </Nav>
      </Header>
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
