import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

export default function Download({ exposureDownload }) {
  const [downloading, setDownload] = useState(false);

  async function downloadPublic() {
    setDownload(true);
    await exposureDownload();
    setDownload(false);
  }

  return (
    <div className="p-4">
      <p>
        Use the following links to download all matrices of mutational signature
        activities from selected study:{' '}
      </p>
      <Button className="p-0" variant="link" onClick={() => downloadPublic()}>
        <LoadingOverlay active={downloading} />
        Download matrixes for selected exposure of mutational signatures and
        study
      </Button>
    </div>
  );
}
