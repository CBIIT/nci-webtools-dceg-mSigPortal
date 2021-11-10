import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../services/store/modal';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

export default function Download({ exposureDownload }) {
  const [downloading, setDownload] = useState(false);
  const [downloadingWorkspace, setWorkspace] = useState(false);

  const dispatch = useDispatch();
  const exposure = useSelector((state) => state.exposure);
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  async function downloadPublic() {
    setDownload(true);
    await exposureDownload();
    setDownload(false);
  }

  // downloads current sesssion from tmp/[id] along with redux state
  async function downloadWorkspace() {
    setWorkspace(true);

    const response = await fetch(`api/downloadWorkspace`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        id: exposure.exposureState.projectID,
        state: {
          ...exposure,
          exposureState: { ...exposure.exposureState, displayTab: 'tmb' },
        },
      }),
    });

    if (response.ok) {
      const objectURL = URL.createObjectURL(await response.blob());
      const tempLink = document.createElement('a');

      tempLink.href = `${objectURL}`;
      tempLink.setAttribute('download', `msigportal-workspace.tgz`);
      document.body.appendChild(tempLink);
      tempLink.click();
      document.body.removeChild(tempLink);
    } else {
      mergeError(`error`);
    }
    setWorkspace(false);
  }

  return (
    <div className="p-4">
      <p>
        Use the following links to download all matrices of mutational signature
        activities from selected study:{' '}
      </p>
      <div>
        <Button variant="link" onClick={() => downloadPublic()}>
          <LoadingOverlay active={downloading} />
          Download matrixes for selected exposure of mutational signatures and
          study
        </Button>
      </div>
      <div>
        <Button variant="link" onClick={() => downloadWorkspace()}>
          <LoadingOverlay active={downloadingWorkspace} />
          Download Workspace
        </Button>
      </div>
    </div>
  );
}
