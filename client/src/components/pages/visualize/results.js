import React, { useState, useEffect } from 'react';
import { Form } from 'react-bootstrap';

const { Group, Label, Control } = Form;

const root =
  process.env.NODE_ENV === 'development'
    ? 'http://localhost:8330/'
    : window.location.pathname;

export default function Results({ plots }) {
  const [svg, setSvg] = useState('');
  const plotDir = plots.plotDir;

  useEffect(() => {
    if (plots.plots) handlePlot(plots.plots[0]);
  }, [plots]);

  const handlePlot = async (plot) => {
    const response = await fetch(`${root}svg`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ path: `${plotDir}/${plot}` }),
    });
    const pic = await response.blob();
    let objectURL = URL.createObjectURL(pic);
    console.log(pic, objectURL);
    setSvg(objectURL);
  };
  return plots.plots ? (
    <div>
      <Form>
        <Group>
          <Label>Choose SVG</Label>
          <Control
            as="select"
            // value={inputFormat}
            onChange={(e) => handlePlot(e.target.value)}
          >
            {plots.plots.map((plot) => {
              return <option value={plot}>{plot}</option>;
            })}
          </Control>
        </Group>
      </Form>
      <img className="w-100" src={svg}></img>
    </div>
  ) : (
    <div>
      <h2>Instructions</h2>
      <p>Upload Sample Variants and specify parameters</p>
    </div>
  );
}
