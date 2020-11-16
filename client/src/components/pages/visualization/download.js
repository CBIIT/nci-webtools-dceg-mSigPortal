import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

export default function Download() {
  const { projectID, downloads, statistics } = useSelector(
    (state) => state.visualizeResults
  );

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
      dispatchError(`${file} is not available`);
    }
    setDownload((downloading) => downloading.filter((item) => item != file));
  }
  return (
    <div>
      {statistics.length > 0 ? (
        <p>{statistics}</p>
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
