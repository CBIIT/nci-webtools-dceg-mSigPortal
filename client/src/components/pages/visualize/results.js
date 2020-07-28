import React, { useEffect } from 'react';
import { Form, Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualizeResults,
} from '../../../services/store';

import PyTab from './pyTab';
import RTab from './rTab';

const { Group, Label, Control } = Form;
const { Header, Body } = Card;
const { Item, Link } = Nav;

const root =
  process.env.NODE_ENV === 'development'
    ? 'http://localhost:8330/'
    : window.location.pathname;

export default function Results() {
  const {
    error,
    projectID,
    displayTab,
    pyPlotURL,
    rPlotURL,
    pyTab,
    rTab,
  } = useSelector((state) => state.visualizeResults);

  const { mapping, filtered } = pyTab;

  const {
    profileType,
    signatureSet,
    selectName2,
    selectSigFormula,
    sigFormula,
    rPlots,
  } = rTab;

  // get mapping after retrieving projectID
  useEffect(() => {
    if (projectID.length) {
      getSummary(projectID);
    }
  }, [projectID]);

  async function setPlot(index, type = 'python') {
    let plot = filtered[index];
    if (type == 'r') plot = rPlots[index];
    if (plot) {
      try {
        const response = await fetch(`${root}visualize/svg`, {
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

          if (pyPlotURL.length) URL.revokeObjectURL(pyPlotURL);
          if (type == 'python') {
            dispatchVisualizeResults({
              pyPlotURL: objectURL,
            });
          } else if (type == 'r') {
            dispatchVisualizeResults({
              rPlotIndex: index,
              rPlotURL: objectURL,
            });
          }
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      console.log('invalid index', index);
    }
  }

  // retrieve mapping of samples to plots from summary file
  async function getSummary() {
    const response = await fetch(`${root}visualize/summary`, {
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
      rTab: { ...rTab, selectName2: nameOptions[0] },
    });
  }

  async function submitR() {
    dispatchVisualizeResults({ rTab: { ...rTab, submitOverlay: true } });
    let args = {
      profileName: profileType,
      signatureSetName: signatureSet,
      sampleName: selectName2,
      signatureName: null,
      formula: null,
      projectID: projectID,
    };
    selectSigFormula == 'signature'
      ? (args.signatureName = sigFormula)
      : (args.formula = sigFormula);

    try {
      const response = await fetch(`${root}api/visualizeR`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(args),
      });
      if (!response.ok) {
        const err = await response.text();
        console.log(err);
        dispatchVisualizeResults({
          rTab: { ...rTab, debugR: err, rPlots: [] },
        });
      } else {
        const data = await response.json();
        console.log(data);
        dispatchVisualizeResults({
          rTab: { ...rTab, debugR: data.debugR, rPlots: data.plots },
        });
      }
    } catch (err) {
      dispatchError(err);
    } finally {
      dispatchVisualizeResults({ rTab: { ...rTab, submitOverlay: false } });
    }
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
              active={displayTab == 'r'}
              onClick={() => dispatchVisualizeResults({ displayTab: 'r' })}
            >
              R
            </Link>
          </Item>
        </Nav>
      </Header>
      <Body style={{ display: displayTab == 'python' ? 'block' : 'none' }}>
        <PyTab setPlot={(e) => setPlot(e)} />
      </Body>
      <Body style={{ display: displayTab == 'r' ? 'block' : 'none' }}>
        <RTab
          setPlot={(e, type) => setPlot(e, type)}
          submitR={() => submitR()}
        />
      </Body>
    </Card>
  ) : (
    <div>
      <h2>Instructions</h2>
      <p>Upload Sample Variants and specify parameters</p>
    </div>
  );
}
