import React, { useEffect } from 'react';
import { Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualizeResults,
  dispatchPyTab,
  dispatchCosineSimilarity,
  dispatchProfileComparison,
  dispatchPCA,
} from '../../../services/store';

import PyTab from './pyTab';
import CosineSimilarity from './cosineSimilarity';
import ProfileComparison from './profileComparison';
import PCA from './pca';
import Download from './download';

const { Header, Body } = Card;
const { Item, Link } = Nav;

export default function Results() {
  const { error, projectID, displayTab, downloads, summary } = useSelector(
    (state) => state.visualizeResults
  );
  const pyTab = useSelector((state) => state.pyTab);
  const rootURL = window.location.pathname;

  // get mapping of plots after retrieving projectID
  useEffect(() => {
    if (summary.length) {
      mapSummary();
    }
  }, [summary]);

  // retrieve mapping of samples to plots from summary file
  async function mapSummary() {
    const nameOptions = [...new Set(summary.map((plot) => plot.Sample_Name))];
    const profileOptions = [
      ...new Set(summary.map((plot) => plot.Profile_Type)),
    ];
    const matrixOptions = [...new Set(summary.map((plot) => plot.Matrix))];
    const tagOptions = [...new Set(summary.map((plot) => plot.Tag))];

    const selectName = pyTab.selectName || nameOptions[0];
    const selectProfile = pyTab.selectProfile || profileOptions[0];
    const selectMatrix = pyTab.selectMatrix || matrixOptions[0];
    const selectTag = pyTab.selectTag || tagOptions[0];

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

    dispatchPyTab({
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
      withinMatrixSize: filteredMatrixOptions[0],
      withinMatrixOptions: filteredMatrixOptions,
    });

    dispatchProfileComparison({
      withinProfileType: profileOptions[0],
      withinSampleName1: nameOptions[0],
      withinSampleName2: nameOptions[1],
      refProfileType: profileOptions[0],
      refSampleName: nameOptions[0],
    });

    dispatchPCA({ profileType: profileOptions[0] });
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
  ) : summary.length ? (
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
      <Body style={{ display: displayTab == 'python' ? 'block' : 'none' }}>
        <PyTab />
      </Body>
      <Body
        style={{ display: displayTab == 'cosineSimilarity' ? 'block' : 'none' }}
      >
        <CosineSimilarity
          downloadResults={(path) => downloadResults(path)}
          submitR={(fn, args) => submitR(fn, args)}
        />
      </Body>
      <Body
        style={{
          display: displayTab == 'profileComparison' ? 'block' : 'none',
        }}
      >
        <ProfileComparison submitR={(fn, args) => submitR(fn, args)} />
      </Body>
      <Body style={{ display: displayTab == 'pca' ? 'block' : 'none' }}>
        <PCA
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
