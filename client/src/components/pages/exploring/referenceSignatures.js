import React, { useEffect } from 'react';
import { Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExploringRefSigs,
} from '../../../services/store';

import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

export default function ReferenceSignatures({ submitR }) {
  const rootURL = window.location.pathname;
  const { plotPath, plotURL, debugR, err, displayDebug, loading } = useSelector(
    (state) => state.exploringRefSigs
  );
  const { displayTab } = useSelector((state) => state.exploring);

  useEffect(() => {
    if (!loading && !plotPath && displayTab == 'referenceSignatures') {
      calculateR('referenceSignatures', {});
    }
  }, [plotPath, displayTab]);

  async function calculateR(fn, args) {
    console.log(fn);
    dispatchExploringRefSigs({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchExploringRefSigs({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output } = await response.json();

        dispatchExploringRefSigs({
          debugR: debugR,
          loading: false,
          plotPath: output.plotPath,
        });
        setRPlot(output.plotPath, 'within');
      }
    } catch (err) {
      dispatchError(err);
      dispatchExploringRefSigs({ loading: false });
    }
  }

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`${rootURL}getSVG`, {
          method: 'POST',
          headers: {
            Accept: 'image/svg',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ path: plotPath }),
        });
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL) URL.revokeObjectURL(plotURL);
          dispatchExploringRefSigs({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExploringRefSigs({ err: true, plotURL: '' });
    }
  }

  return (
    <div>
      <LoadingOverlay active={loading} />
      <div id="withinPlot">
        <div style={{ display: err ? 'block' : 'none' }}>
          <p>An error has occured. Check the debug section for more info.</p>
        </div>
        <div style={{ display: plotURL ? 'block' : 'none' }}>
          <div className="d-flex">
            <a
              className="px-2 py-1"
              href={plotURL}
              download={plotURL.split('/').slice(-1)[0]}
            >
              Download Plot
            </a>
          </div>
          <div className="p-2 border rounded">
            <Row>
              <Col>
                <img className="w-100 my-4 h-1000" src={plotURL}></img>
              </Col>
            </Row>
          </div>
        </div>
      </div>

      <Button
        variant="link"
        className="p-0 mt-5"
        onClick={() =>
          dispatchExploringRefSigs({
            displayDebug: !displayDebug,
          })
        }
      >
        R Debug
      </Button>
      <pre
        className="border rounded p-1 "
        style={{ display: displayDebug ? 'block' : 'none' }}
      >
        <div className="border">
          {Array.isArray(debugR) ? (
            debugR.map((line, index) => {
              return (
                <p key={index} className="m-0">
                  [{index}] {line}
                </p>
              );
            })
          ) : (
            <p>{debugR}</p>
          )}
        </div>
      </pre>
    </div>
  );
}
