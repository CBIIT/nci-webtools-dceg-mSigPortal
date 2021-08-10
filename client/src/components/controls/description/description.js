import React, { useState } from 'react';
import { Button } from 'react-bootstrap';

export default function Debug({ short, description, ...props }) {
  const [show, setShow] = useState(false);

  return (
    <p {...props}>
      {!show ? short : description}{' '}
      <span>
        <Button
          className="p-0 ml-3"
          variant="link"
          onClick={() => setShow(!show)}
        >
          Show {show ? 'Less' : 'More'}
        </Button>
      </span>
    </p>
  );
}
