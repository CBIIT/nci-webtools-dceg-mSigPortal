import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
import { saveAs } from 'file-saver';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../services/store/modal';

export default function Download() {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    source,
    study,
    experimentalStrategy,
    cancerType,
    id,
    downloads,
    statistics,
  } = visualization.main;

  const [downloading, setDownload] = useState([]);
  const [downloadingWorkspace, setWorkspace] = useState(false);

  async function downloadOutput(file) {
    setDownload((downloading) => [...downloading, file]);
    const response = await fetch(
      `web/visualization/download?id=${id}&file=${file}`
    );
    if (response.ok) {
      saveAs(
        await response.blob(),
        file.split('/')[file.split('/').length - 1]
      );
    } else {
      mergeError(`${file} is not available`);
    }
    setDownload((downloading) => downloading.filter((item) => item != file));
  }

  async function downloadPublic() {
    setDownload([1]);
    const response = await fetch(`web/visualization/downloadPublic`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        fn: 'downloadPublicData',
        args: {
          id,
          study: study,
          cancerType: cancerType,
          experimentalStrategy: experimentalStrategy,
        },
      }),
    });

    if (response.ok) {
      saveAs(
        await response.blob(),
        `msigportal-${study}-${cancerType}-${experimentalStrategy}.tar.gz`
      );
    } else {
      mergeError(`public data is not available`);
    }
    setDownload([]);
  }

  // downloads current sesssion from tmp/[id] along with redux state
  async function downloadWorkspace() {
    setWorkspace(true);

    const { mutationalProfiles, ...rest } = visualization;
    const response = await fetch(`web/downloadWorkspace`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        id,
        state: {
          visualization: {
            ...rest,
            state: {
              ...rest.state,
              email: '',
              queueMode: false,
              displayTab: 'profilerSummary',
            },
          },
        },
      }),
    });

    if (response.ok) {
      saveAs(await response.blob(), 'visualization-workspace.zip');
    } else {
      mergeError(
        `error - Please calculate Cosine Similarity tab before attempt to download the workspace`
      );
    }
    setWorkspace(false);
  }

  return (
    <div className="bg-white border rounded p-4">
      {source == 'user' ? (
        <>
          {statistics.length > 0 ? (
            <div>
              <p>{statistics}</p>
              <p>
                Use the following links to download all analyzed data and
                visualization from selected or input dataset:
              </p>
            </div>
          ) : (
            <p>No statistics available</p>
          )}
          {downloads.length > 0 ? (
            <div>
              {downloads.map((file) => (
                <div key={file}>
                  <Button variant="link" onClick={() => downloadOutput(file)}>
                    <LoadingOverlay active={downloading.indexOf(file) != -1} />
                    Download {file.split('/')[file.split('/').length - 1]}
                  </Button>
                </div>
              ))}
            </div>
          ) : (
            <p>No files available</p>
          )}
        </>
      ) : (
        <div>
          <Button variant="link" onClick={() => downloadPublic()}>
            <LoadingOverlay active={downloading.length > 0} />
            Download matrices of different mutational profiles for selected
            study
          </Button>
        </div>
      )}
      <Button variant="link" onClick={() => downloadWorkspace()}>
        <LoadingOverlay active={downloadingWorkspace} />
        Download Workspace
      </Button>
    </div>
  );
}
