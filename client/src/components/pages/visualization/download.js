import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
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
    projectID,
    downloads,
    statistics,
  } = visualization.state;

  const [downloading, setDownload] = useState([]);
  const [downloadingSession, setSession] = useState(false);

  async function downloadOutput(file) {
    setDownload((downloading) => [...downloading, file]);
    const response = await fetch(
      `api/visualization/download?id=${projectID}&file=${file}`
    );
    if (response.ok) {
      const objectURL = URL.createObjectURL(await response.blob());
      const tempLink = document.createElement('a');

      tempLink.href = `${objectURL}`;
      tempLink.setAttribute('download', `${file}`);
      document.body.appendChild(tempLink);
      tempLink.click();
      document.body.removeChild(tempLink);
    } else {
      mergeError(`${file} is not available`);
    }
    setDownload((downloading) => downloading.filter((item) => item != file));
  }

  async function downloadPublic() {
    setDownload([1]);
    const response = await fetch(`api/visualization/downloadPublic`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        fn: 'downloadPublicData',
        args: {
          id: projectID,
          study: study,
          cancerType: cancerType,
          experimentalStrategy: experimentalStrategy,
        },
      }),
    });

    if (response.ok) {
      const objectURL = URL.createObjectURL(await response.blob());
      const tempLink = document.createElement('a');

      tempLink.href = `${objectURL}`;
      tempLink.setAttribute(
        'download',
        `msigportal-${study}-${cancerType}-${experimentalStrategy}`
      );
      document.body.appendChild(tempLink);
      tempLink.click();
      document.body.removeChild(tempLink);
    } else {
      mergeError(`public data is not available`);
    }
    setDownload([]);
  }

  // downloads current sesssion from tmp/[id] along with redux state
  async function downloadSession() {
    setSession(true);

    const response = await fetch(`api/downloadSession`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        id: projectID,
        state: {
          visualization: {
            ...visualization,
            state: {
              ...visualization.state,
              email: '',
              queueMode: false,
              displayTab: 'profilerSummary',
            },
            mutationalProfiles: {
              ...visualization.mutationalProfiles,
              plotPath: '',
            },
          },
        },
      }),
    });

    if (response.ok) {
      const objectURL = URL.createObjectURL(await response.blob());
      const tempLink = document.createElement('a');

      tempLink.href = `${objectURL}`;
      tempLink.setAttribute('download', `msigportal-session.tgz`);
      document.body.appendChild(tempLink);
      tempLink.click();
      document.body.removeChild(tempLink);
    } else {
      mergeError(`error`);
    }
    setSession(false);
  }

  return (
    <div className="bg-white border rounded p-4">
      {source == 'user' ? (
        <>
          {statistics.length > 0 ? (
            <div>
              <p>{statistics}</p>
              <p>
                Using the following links to download different mutational
                profiles, plots and the mapping information between mutation and
                each mutation type.
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
                    Download {file.split('.')[0]}
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
      <Button variant="link" onClick={() => downloadSession()}>
        <LoadingOverlay active={downloadingSession} />
        Download Session
      </Button>
    </div>
  );
}
