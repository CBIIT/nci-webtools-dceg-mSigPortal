import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';

const { Label, Group } = Form;

export default function MsLandscape({ calculateLandscape, handleVariable }) {
  const exposure = useSelector((state) => state.exposure);
  const { variableFile, plotPath, debugR, err, loading } = exposure.msLandscape;
  const { projectID, source } = exposure.exposureState;

  return (
    <div>
      <div className="p-3">
        <b>Landscape of Mutational Signature Activity</b>
        <p className="m-0">
          This page allows for you investigate the overall landscape of exposure
          of mutational signatures across samples of a given cancer type. In the
          combined plot below, you will find a series of bar plots on top of one
          another to form the overall plot. The x-axis is the samples for the
          cancer type and study selected from the left panel. Starting with the
          bottom plot, this illustrates the contribution of different signatures
          within a sample. The colors that correspond to each mutational
          signature can be found at the bottom in the “Mutational Signatures”
          legend. Just above the contributions plot is a cosine similarity bar
          of the mutation signatures deconvolution found in each sample. The
          plot above the cosine similarity assignments is a plot that contains
          the number of mutations in each sample assigned to each mutational
          signature. The topmost plot is unsupervised clustering of the samples
          based on signature contributions. Additional bars can add under the
          clustering by uploading the variable data.
        </p>
      </div>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col lg="4">
            <Group controlId="landscape">
              <Label>Upload Variable Data</Label>
              <Form.File
                id="variableData"
                label={variableFile || 'Upload here (optional)'}
                // accept=''
                onChange={(e) => handleVariable(e.target.files[0])}
                custom
              />
            </Group>
          </Col>
          <Col lg="1" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="secondary"
              onClick={() => {
                handleVariable(new File([], ''));
              }}
              disabled={!variableFile}
            >
              Remove
            </Button>
          </Col>
          <Col />
          <Col lg="2" className="d-flex">
            <Button
              disabled={source == 'user' && !projectID}
              className="ml-auto mb-auto"
              variant="primary"
              onClick={calculateLandscape}
            >
              Recalculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="exposureLandscapePlot">
        {err && (
          <div>
            <hr />
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        {plotPath && (
          <>
            <hr />
            <Plot
              className="p-3"
              title="Landscape of Mutational Signature Activity"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${projectID}${plotPath}`}
              height="1500px"
            />
          </>
        )}
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
