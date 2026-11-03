import { useForm } from 'react-hook-form';
import MsLandscapeForm from './msLandscape-form';
import MsLandscapePlot from './msLandscape-plot';

// one dataset's variable-file upload + plot, so Primary/Comparison each keep their own upload
export default function MsLandscapePanel({ state, fieldId, divId }) {
  const { control, watch, setValue, resetField } = useForm({
    defaultValues: { variableFile: '' },
  });
  const { variableFile } = watch();

  return (
    <>
      <MsLandscapeForm
        state={state}
        variableFile={variableFile}
        formControl={control}
        setValue={setValue}
        resetField={resetField}
        fieldId={fieldId}
      />
      <hr />
      <MsLandscapePlot state={state} variableFile={variableFile} divId={divId} />
    </>
  );
}
