import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/modal';
import { getBlob } from '../../../../services/utils';
import { saveAs } from 'file-saver';

export default function Download() {
  const dispatch = useDispatch();

  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const [downloading, setDownload] = useState(false);

  async function downloadPublic() {
    setDownload(true);
    try {
      saveAs(
        await getBlob('Signature/signature_refsets.txt.zip'),
        'signature_refsets.txt.zip'
      );
    } catch (err) {
      console.log(err);
      mergeError(`public data is not available`);
    }
    setDownload(false);
  }

  return (
    <div className="p-4">
      <p>
        Use the following link to download all matrices of collected mutational
        signatures:
      </p>
      <Button className="p-0" variant="link" onClick={() => downloadPublic()}>
        <LoadingOverlay active={downloading.length > 0} />
        Download matrixes for all available mutational signatures in mSigPortal
      </Button>
    </div>
  );
}
