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
  };

  return (
    <AsyncCreatableSelect
      name={name}
      placeholder={placeholder}
      value={value}
      onChange={onChange}
      isMulti
      formatCreateLabel={(userInput) => `${placeholder}: ${userInput}`}
      {...customStyle}
      {...rest}
    />
  );
}
