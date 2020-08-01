import React, { useEffect } from 'react';
import { Form, Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualizeResults,
} from '../../../services/store';

import PyTab from './pyTab';
import { CosineSimilarity } from './rTab';

const { Header, Body } = Card;
const { Item, Link } = Nav;

export default function Results() {
  const { error, projectID, displayTab, pyTab, cosineSimilarity } = useSelector(
    (state) => state.visualizeResults
  );

  const { mapping } = pyTab;

  // get mapping after retrieving projectID
  useEffect(() => {
    if (projectID.length) {
      getSummary(projectID);
    }
  }, [projectID]);

  // retrieve mapping of samples to plots from summary file
  async function getSummary() {
    const response = await fetch(`/visualize/summary`, {
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

    const selectName = pyTab.selectName || nameOptions[0];
    const selectProfile = pyTab.selectProfile || profileOptions[0];
    const selectMatrix = pyTab.selectMatrix || matrixOptions[0];
    const selectTag = pyTab.selectTag || tagOptions[0];

    const filteredPlots = mapping.filter(
      (plot) =>
        plot.Sample_Name == selectName &&
        plot.Profile_Type == selectProfile &&
        plot.Matrix == selectMatrix &&
        plot.Tag == selectTag
    );

    const filteredProfileOptions = [
      ...new Set(
        mapping
          .filter((plot) => plot.Sample_Name == selectName)
          .map((plot) => plot.Profile_Type)
      ),
    ];
    const filteredMatrixOptions = [
      ...new Set(
        mapping
          .filter((plot) => plot.Profile_Type == selectProfile)
          .map((plot) => plot.Matrix)
      ),
    ];
    const filteredTagOptions = [
      ...new Set(
        mapping
          .filter((plot) => plot.Matrix == selectMatrix)
          .map((plot) => plot.Tag)
      ),
    ];

    dispatchVisualizeResults({
      pyTab: {
        ...pyTab,
        mapping: mapping,
        filtered: filteredPlots,
        nameOptions: nameOptions,
        profileOptions: filteredProfileOptions,
        matrixOptions: filteredMatrixOptions,
        tagOptions: filteredTagOptions,
        selectName: selectName,
        selectProfile: selectProfile,
        selectMatrix: selectMatrix,
        selectTag: selectTag,
      },
      cosineSimilarity: {
        ...cosineSimilarity,
        profileType1: profileOptions[0],
        profileType2: profileOptions[0],
        matrixSize: matrixOptions[0],
        matrixOptions: matrixOptions,
      },
    });
  }

  function submitR(fn, args) {
    return fetch(`/api/visualizeR`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ fn: fn, args: args, projectID: projectID }),
    });
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
              active={displayTab == 'cosineSimilarity'}
              onClick={() =>
                dispatchVisualizeResults({ displayTab: 'cosineSimilarity' })
              }
            >
              Cosine Similarity
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
        <CosineSimilarity submitR={(fn, args) => submitR(fn, args)} />
      </Body>
    </Card>
  ) : (
    <div>
      <h2>Instructions</h2>
      <p>Upload Sample Variants and specify parameters</p>
    </div>
  );
}
