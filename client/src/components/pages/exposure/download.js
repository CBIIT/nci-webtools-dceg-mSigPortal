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
      <Button variant="link" onClick={() => downloadPublic()}>
        <LoadingOverlay active={downloading} />
        Download matrixes for selected exposure of mutational signatures and
        study
      </Button>
    </div>
  );
}
