import { Container, Nav, Tab } from 'react-bootstrap';
import MsLandscapePanel from './msLandscape-panel';
import Description from '@/components/controls/description/description';
import { useState } from 'react';
import { getDatasetLabel } from '@/services/utils';

export default function MsLandscape({ state, compareState }) {
  const [display, setDisplay] = useState('primary');

  return (
    <Container fluid className="bg-white border rounded p-0">
      <div className="p-3">
        <b>Landscape of Mutational Signature Activity</b>
        <Description
          className="m-0"
          less="The following clustering of mutational signatures allows users to investigate the overall landscape of mutational signatures from the selected cancer type."
          more="The combined plot below includes the following sub-plots from top to bottom: 1) Unsupervised clustering of all samples based on signature contributions; 2) Values from uploaded variable data assigned to each sample (optional); 3) Stacked bar plot that contains the number of mutations assigned to each mutational signature in each sample; 4) Cosine similarity between the original mutational pattern and reconstructed mutational pattern for each sample, which indicates the performance of mutational signature deconvolution; 5) Signature contribution plot, which illustrates the contribution of different signatures within a sample. The colors that correspond to each subplot can be found in the legend."
        />
      </div>
      <hr />
      {compareState ? (
        <Tab.Container
          transition={false}
          activeKey={display}
          onSelect={(tab) => setDisplay(tab)}
        >
          <Nav variant="tabs" className="px-3 mt-3">
            <Nav.Item>
              <Nav.Link
                eventKey="primary"
                as="button"
                className="outline-none exploration"
              >
                <strong>{getDatasetLabel(state) || 'Primary Data'}</strong>
              </Nav.Link>
            </Nav.Item>
            <Nav.Item>
              <Nav.Link
                eventKey="compare"
                as="button"
                className="outline-none exploration"
              >
                <strong>{getDatasetLabel(compareState)}</strong>
              </Nav.Link>
            </Nav.Item>
          </Nav>
          <Tab.Content
            className="bg-white tab-pane-bordered rounded-0 d-block mx-3 mb-3"
            style={{ overflowX: 'auto' }}
          >
            <Tab.Pane eventKey="primary" mountOnEnter className="border-0">
              <MsLandscapePanel
                state={state}
                fieldId="variableFile"
                divId="msLandscape-primary"
              />
            </Tab.Pane>
            <Tab.Pane eventKey="compare" mountOnEnter className="border-0">
              <MsLandscapePanel
                state={compareState}
                fieldId="compareVariableFile"
                divId="msLandscape-compare"
              />
            </Tab.Pane>
          </Tab.Content>
        </Tab.Container>
      ) : (
        <MsLandscapePanel
          state={state}
          fieldId="variableFile"
          divId="msLandscape-primary"
        />
      )}
    </Container>
  );
}
