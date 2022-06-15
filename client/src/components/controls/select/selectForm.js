import React from 'react';
import { Form } from 'react-bootstrap';
import ReactSelect, { createFilter } from 'react-select';

export default function Select({
  className,
  id,
  label,
  value,
  options,
  onChange,
  disabled,
  labelClass,
  ...rest
}) {
  const selectStyles = {
    styles: {
      menuPortal: (base) => ({ ...base, zIndex: 9999 }),
      container: (base) => ({
        ...base,
        flex: 1,
      }),
      singleValue: ({
        maxWidth,
        position,
        top,
        transform,
        ...otherStyles
      }) => ({ ...otherStyles }),
      menu: (base) => ({
        ...base,
        width: 'max-content',
        minWidth: '100%',
      }),
    },
    menuPortalTarget: document.body,
    filterOption: createFilter({ ignoreAccents: false }),
  };

  return (
    <Form.Group controlId={id} className={className}>
      {label && <Form.Label className={labelClass}>{label}</Form.Label>}
      <ReactSelect
        name={id}
        inputId={id}
        options={options}
        value={value}
        onChange={onChange}
        isDisabled={disabled}
        {...selectStyles}
        {...rest}
      />
    </Form.Group>
  );
}
