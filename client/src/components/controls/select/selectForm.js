import React from 'react';
import { Form } from 'react-bootstrap';
import ReactSelect, { createFilter } from 'react-select';
import { Controller } from 'react-hook-form';

export default function Select({
  name,
  label,
  options,
  disabled,
  className,
  labelClass,
  control,
  rules,
  defaultValue,
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
    <Form.Group controlId={name} className={className}>
      {label && <Form.Label className={labelClass}>{label}</Form.Label>}
      <Controller
        name={name}
        control={control}
        rules={rules}
        defaultValue={defaultValue}
        render={({ field }) => (
          <ReactSelect
            {...selectStyles}
            {...field}
            name={name}
            inputId={name}
            options={options}
            isDisabled={disabled}
            {...rest}
          />
        )}
      />
    </Form.Group>
  );
}
