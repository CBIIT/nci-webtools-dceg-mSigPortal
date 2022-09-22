import Description from '../../../../controls/description/description';
import RsProfileFormPlot from '../rsProfile/rsProfile-form-plot';
import { useSelector, useDispatch } from 'react-redux';

export default function RsProfile() {
  const store = useSelector((state) => state.catalog);
  const { plots } = store.rSProfiles;
  return (
    <div>
      <Description
        className="p-3 m-0"
        less="Enter any [Signature Source], [Profile Name], [Reference Signature Set], [Experimental Strategy], and [Signature Name] below to visualize the mutational signature profile."
        more="Click ‘+ Add Plot’ to load one or more mutational signature profiles at the same time."
      />
      <hr />
      {plots.map((e, i) => (
        <RsProfileFormPlot options={e} index={i} />
      ))}
    </div>
  );
}