import React, { useEffect } from 'react';
import { Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchProfilerSummary,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

export default function ProfilerSummary({ submitR }) {
  const {
    source,
    study,
    studyOptions,
    cancerType,
    pubExperimentalStrategy,
    pDataOptions,
  } = useSelector((state) => state.visualize);
  const { matrixList, svgList } = useSelector(
    (state) => state.visualizeResults
  );
  const { plotPath, plotURL, err, debugR, displayDebug, loading } = useSelector(
    (state) => state.profilerSummary
  );
  const rootURL = window.location.pathname;

  useEffect(() => {
    if (!loading && !plotPath) {
      if (source == 'user' && matrixList.length) {
        calculateR('profilerSummary', {
          matrixList: JSON.stringify(matrixList),
        });
      } else if (source == 'public' && svgList.length) {
        calculateR('profilerSummaryPublic', {
          study: study,
          cancerType: cancerType,
          experimentalStrategy: pubExperimentalStrategy,
        });
      }
    }
  }, [matrixList, svgList]);

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
          dispatchProfilerSummary({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchProfilerSummary({ err: true, plotURL: '' });
    }
  }

  async function calculateR(fn, args) {
    console.log(fn);
    dispatchProfilerSummary({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchProfilerSummary({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output } = await response.json();

        dispatchProfilerSummary({
          debugR: debugR,
          loading: false,
          plotPath: output.plotPath,
        });
        setRPlot(output.plotPath, 'within');
      }
    } catch (err) {
      dispatchError(err);
      dispatchProfilerSummary({ loading: false });
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
              download={plotPath.split('/').slice(-1)[0]}
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
          dispatchProfilerSummary({
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
