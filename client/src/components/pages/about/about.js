import React from 'react';
import { useDispatch, useSelector } from 'react-redux';
import { updateVisualize } from '../../../services/actions';
export default function About() {
  const dispatch = useDispatch();
  const {
    count
  } = useSelector(state => state.visualize);
  
  return (
    <div
      className="mt-3 container bg-white tab-pane-bordered rounded-0 p-4"
      style={{
        minHeight: '420px',
      }}
    >
      <h1 className="font-weight-light">About mSigPortal</h1>
      <hr />

      <p>TBA</p>

      <button
        onClick={e => {
          dispatch(updateVisualize({ count: count + 1 }));
        }}>
        INCREMENT
      </button>
      <button 
        onClick={e => {
          dispatch(updateVisualize({ count: count - 1}));
        }}>
        DECREMENT
      </button>
      <pre>
        { count }
      </pre>

    </div>
  );
}
