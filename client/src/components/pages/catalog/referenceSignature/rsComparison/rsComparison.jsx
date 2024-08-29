import Description from '../../../../controls/description/description';
import RsComparisonPlot from './form-plot';

export default function RsComparison() {
  return (
    <div>
      <Description
        className="p-3 m-0"
        less="Below you can compare two mutational signatures from curated reference signature sets. "
        more="Use the dropdown menus to input a [Profile Name], two [Reference Signature Sets], and two [Signature Names] within the selected Signature Sets."
      />
      <hr />
      <RsComparisonPlot />
    </div>
  );
}
