import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
import { saveAs } from 'file-saver';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../services/store/modal';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

export default function Download({ exposureDownload }) {
  const [downloading, setDownload] = useState(false);
  const [downloadingWorkspace, setWorkspace] = useState(false);

  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
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

    const response = await fetch(`web/downloadWorkspace`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        id: exploration.main.id,
        state: {
          ...exploration,
          main: { ...exploration.main, displayTab: 'tmb' },
        },
      }),
    });

    if (response.ok) {
      saveAs(await response.blob(), 'exploration-workspace.zip');
    } else {
      mergeError(`error`);
    }
    setWorkspace(false);
  }

  return (
    <div className="p-4">
      <p>
        Use the following links to download all matrices of mutational signature
        activities from the selected study:{' '}
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
