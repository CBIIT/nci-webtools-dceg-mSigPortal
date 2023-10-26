import React from 'react';
import { Form } from 'react-bootstrap';
import ReactSelect, { createFilter } from 'react-select';

const { Group, Label } = Form;

export default function CustomSelect({
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
  const props = {
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
    getOptionLabel: ({ label }) => (label == 'NA' ? 'N/A' : label),
    filterOption: createFilter({ ignoreAccents: false }),
  };

  // parse array of strings into array of option objects
  // non array args return an empty array
  const optionsObject = Array.isArray(options)
    ? typeof options[0] != 'object'
      ? options.map((v) => ({ value: v, label: v }))
      : options
    : [];

  return (
    <Group controlId={id} className={className}>
      {label && <Label className={labelClass}>{label}</Label>}
      <ReactSelect
        name={id}
        inputId={id}
        options={optionsObject}
        value={optionsObject.filter((option) => option.value === value)}
        onChange={(option) => onChange(option.value)}
        isDisabled={disabled}
        {...props}
        {...rest}
      />
    </Group>
  );
}
