import React from 'react';
import { Form } from 'react-bootstrap';
import ReactSelect, { createFilter } from 'react-select';
import { Controller } from 'react-hook-form';

export default function SelectForm({
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
        // override border radius to match the box
        borderRadius: 0,
        // kill the gap
        marginTop: 0,
      }),
      control: (base, state) => ({
        ...base,
        background: '#f1e4ef',
        // match with the menu
        borderRadius: state.isFocused ? '3px 3px 0 0' : 3,
        // Overwrittes the different states of border
        borderColor: state.isFocused ? '#f1e4ef' : '#8e4b86',
        // Removes weird border around container
        boxShadow: state.isFocused ? null : null,
        '&:hover': {
          // Overwrittes the different states of border
          borderColor: state.isFocused ? '#8e4b86' : '#f1e4ef',
        },
      }),
      menuList: (base) => ({
        ...base,
        // kill the white space on first and last option
        padding: 0,
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
