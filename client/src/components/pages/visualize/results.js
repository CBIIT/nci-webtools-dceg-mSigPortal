import React, { useState, useEffect } from 'react';
import { Form } from 'react-bootstrap';
import { useDispatch, useSelector } from 'react-redux';
import { updateVisualizeResults } from '../../../services/actions';
import { getInitialState } from '../../../services/store';

const { Group, Label, Control } = Form;

const root =
  process.env.NODE_ENV === 'development'
    ? 'http://localhost:8330/'
    : window.location.pathname;

export default function Results() {
  const dispatch = useDispatch();
  const { projectID, mapping, plots, displayedPlot, plotURL } = useSelector(
    (state) => state.visualizeResults
  );

  useEffect(() => {
    if (projectID.length) {
      getPlotMapping(projectID);
      console.log('useEffect projectID', projectID);
    }
  }, [projectID]);

  useEffect(() => {
    if (plots.length && displayedPlot.length == 0) {
      setPlot(mapping[0].Sample_Name);
      console.log('useEffect mapping');
    }
  }, [mapping]);

  async function getPlotMapping() {
    const response = await fetch(`${root}visualizeSummary`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ projectID: projectID }),
    });
    const mapping = await response.json();

    dispatch(
      updateVisualizeResults({
        mapping: mapping,
        plots: mapping.map((plot) => plot.Sample_Name).sort(),
      })
    );
  }

  async function setPlot(plotName) {
    const plotPath = mapping.filter((plot) => {
      return plot.Sample_Name == plotName;
    });
    const response = await fetch(`${root}svg`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ path: plotPath[0].Location }),
    });
    const pic = await response.blob();
    let objectURL = URL.createObjectURL(pic);
    console.log(pic, objectURL);
    dispatch(
      updateVisualizeResults({
        displayedPlot: plotName,
        plotURL: objectURL,
      })
    );
  }

  return plots.length ? (
    <div>
      <Form>
        <Group>
          <Label>Choose SVG</Label>
          <Control
            as="select"
            value={displayedPlot}
            onChange={(e) => setPlot(e.target.value)}
          >
            {plots.map((plot) => {
              return (
                <option key={plot} value={plot}>
                  {plot}
                </option>
              );
            })}
          </Control>
        </Group>
      </Form>
      <img className="w-100" src={plotURL}></img>
    </div>
  ) : (
    <div>
      <h2>Instructions</h2>
      <p>Upload Sample Variants and specify parameters</p>
    </div>
  );
}
