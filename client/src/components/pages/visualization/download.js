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

  const { projectID, downloads, statistics } = visualization.results;
  const [downloading, setDownload] = useState([]);

  async function downloadOutput(file) {
    setDownload((downloading) => [...downloading, file]);
    const response = await fetch(
      `api/visualize/download?id=${projectID}&file=${file}`
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
  return (
    <div className="bg-white border rounded p-4">
      {statistics.length > 0 ? (
        <div>
          <p>{statistics}</p>
          <p>
            Using the following links to download different mutational profiles,
            plots and the mapping information between mutation and each mutation
            type.
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
    </div>
  );
}
