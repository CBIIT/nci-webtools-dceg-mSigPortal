import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/modal';

export default function Download() {
  const dispatch = useDispatch();

  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const [downloading, setDownload] = useState(false);

  async function downloadPublic() {
    setDownload(true);
    try {
      const response = await fetch(`api/getFileS3`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          path: 'Signature/signature_refsets.txt.tar.gz',
        }),
      });

      if (response.ok) {
        const objectURL = URL.createObjectURL(await response.blob());
        const tempLink = document.createElement('a');

        tempLink.href = `${objectURL}`;
        tempLink.setAttribute('download', `signature_refsets.txt.tar.gz`);
        document.body.appendChild(tempLink);
        tempLink.click();
        document.body.removeChild(tempLink);
      } else {
        mergeError(`public data is not available`);
      }
    } catch (err) {
      console.log(err);
      mergeError(`public data is not available`);
    }
    setDownload(false);
  }

  return (
    <div className="p-4">
      <Button variant="link" onClick={() => downloadPublic()}>
        <LoadingOverlay active={downloading.length > 0} />
        Download matrixes for all available mutational signatures in mSigPortal
      </Button>
    </div>
  );
}
