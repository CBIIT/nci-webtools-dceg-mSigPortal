import React from 'react';
import { Form } from 'react-bootstrap';
import ReactSelect from 'react-select';
const { Group, Label } = Form;

export default function Select({
  id,
  label,
  value,
  options,
  onChange,
  disabled,
}) {
  const selectFix = {
    styles: {
      menuPortal: (base) => ({ ...base, zIndex: 9999 }),
    },
    menuPortalTarget: document.body,
    getOptionLabel: (option) => (option == 'NA' ? 'N/A' : option),
    getOptionValue: (option) => option,
  };

  return (
    <Group controlId={id}>
      <Label>{label}</Label>
      <ReactSelect
        inputId={id}
        options={options}
        value={[value]}
        onChange={onChange}
        isDisabled={disabled}
        {...selectFix}
      />
    </Group>
  );
}
