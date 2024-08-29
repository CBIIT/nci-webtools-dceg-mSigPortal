import React, { useState } from 'react';
import { Button } from 'react-bootstrap';

export default function Debug({ less = '', more = '', ...props }) {
  const [show, setShow] = useState(false);

  return (
    <div {...props}>
      {show ? (
        <>
          {less} {more}
        </>
      ) : (
        less
      )}{' '}
      <span className="d-inline-flex">
        <Button
          className="p-0 border-0"
          variant="link"
          onClick={() => setShow(!show)}
        >
          Show {show ? 'Less' : 'More'}
        </Button>
      </span>
    </div>
  );
}
