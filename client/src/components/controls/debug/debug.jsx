import React, { useState } from 'react';
import { Button } from 'react-bootstrap';

export default function Debug({ msg }) {
  const [display, setDisplay] = useState(false);

  return (
    <div style={{ display: msg.length ? 'block' : 'none' }}>
      <Button
        variant="link"
        className="p-0 mt-3"
        onClick={() => setDisplay(!display)}
      >
        Debug
      </Button>
      <pre
        className="border rounded p-1"
        style={{ display: display ? 'block' : 'none' }}
      >
        <div className="border">
          {Array.isArray(msg) ? (
            msg.map((line, index) => {
              return (
                <p key={index} className="m-0">
                  [{index}] {line}
                </p>
              );
            })
          ) : (
            <p>{msg}</p>
          )}
        </div>
      </pre>
    </div>
  );
}
