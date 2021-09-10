import React from 'react';
import { Form } from 'react-bootstrap';
import ReactSelect, { createFilter } from 'react-select';

const { Group, Label } = Form;

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
    },
    menuPortalTarget: document.body,
    getOptionLabel: (option) =>
      option.label || (option == 'NA' ? 'N/A' : option),
    getOptionValue: (option) => option.value || option,
    filterOption: createFilter({ ignoreAccents: false }),
  };

  return (
    <Group controlId={id} className={className}>
      {label && <Label className={labelClass}>{label}</Label>}
      <ReactSelect
        inputId={id}
        options={options}
        value={[value] || [options[0]]}
        onChange={onChange}
        isDisabled={disabled}
        {...props}
        {...rest}
      />
    </Group>
  );
}
