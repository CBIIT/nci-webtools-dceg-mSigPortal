// import CreatableSelect from 'react-select/creatable';
import AsyncCreatableSelect from 'react-select/async-creatable';

export default function MultiSelect({
  name,
  placeholder,
  value,
  onChange,
  ...rest
}) {
  const customStyle = {
    //   hide prompts and indicators
    noOptionsMessage: () => null,
    components: {
      DropdownIndicator: () => null,
      IndicatorSeparator: () => null,
    },
    styles: {
      control: (provided, state) => ({
        ...provided,
        borderRadius: '2rem',
        borderColor: 'rgb(111, 208, 178)',
        boxShadow: state.isFocused
          ? '0 0 0 0.25rem rgba(111, 208, 178, 0.25)'
          : 'none',
        ':hover': {
          borderColor: 'rgb(111, 208, 178)',
        },
      }),
      multiValue: (provided, state) => ({
        ...provided,
        borderRadius: '.75rem',
      }),
      multiValueRemove: (provided, state) => ({
        ...provided,
        borderRadius: '0 .75rem .75rem 0',
      }),
    },
  };

  return (
    <AsyncCreatableSelect
      name={name}
      placeholder={placeholder}
      value={value}
      onChange={onChange}
      isMulti
      formatCreateLabel={(userInput) => `${placeholder}: ${userInput}`}
    //   {...customStyle}
      {...rest}
    />
  );
}
