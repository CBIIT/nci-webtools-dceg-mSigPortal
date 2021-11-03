import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import Description from '../../controls/description/description';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faFolderMinus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';

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
          <Description
            less="The following clustering of mutational signatures allows users to investigate the overall landscape of mutational signatures from the selected cancer type."
            more="The combined plot below includes the following sub-plots from top to bottom: 1) Unsupervised clustering of all samples based on signature contributions; 2) Values from uploaded variable data assigned to each sample (optional); 3) Stacked bar plot that contains the number of mutations assigned to each mutational signature in each sample; 4) Cosine similarity between the original mutational pattern and reconstructed mutational pattern for each sample, which indicates the performance of mutational signature deconvolution; 5) Signature contribution plot, which illustrates the contribution of different signatures within a sample. The colors that correspond to each subplot can be found in the legend."
          />
        </p>
      </div>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col lg="auto">
            <Group controlId="landscape">
              <Label>Upload Variable Data</Label>
              <div className="d-flex">
                <Form.File
                  id="variableData"
                  label={variableFile || 'Upload here (optional)'}
                  title={variableFile || 'Upload here (optional)'}
                  value={''}
                  // accept=''
                  onChange={(e) => handleVariable(e.target.files[0])}
                  custom
                />
                {variableFile && (
                  <Button
                    className="ml-1"
                    size="sm"
                    title="Remove"
                    variant="danger"
                    disabled={loading}
                    onClick={() => {
                      handleVariable(new File([], ''));
                    }}
                  >
                    <FontAwesomeIcon icon={faFolderMinus} size="lg" />
                  </Button>
                )}
              </div>
            </Group>
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              disabled={source == 'user' && !projectID}
              className="mt-auto mb-3"
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
              plotPath={`api/results/${plotPath}`}
              height="1500px"
            />
          </>
        )}
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
