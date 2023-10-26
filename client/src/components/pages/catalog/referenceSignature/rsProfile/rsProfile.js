import Description from '../../../../controls/description/description';
import RsProfileFormPlot from '../rsProfile/rsProfile-form-plot_1';
import { useSelector, useDispatch } from 'react-redux';

export default function RsProfile() {
  const store = useSelector((state) => state.catalog);
  const { plots } = store.rSProfiles;
  return (
    <div>
      Enter any [Signature Source], [Profile Name], [Reference Signature Set],
      [Experimental Strategy], and [Signature Name] below to visualize the
      mutational signature profile.
      <hr />
      {/* {plots.map((e, i) => (
        <RsProfileFormPlot options={e} index={i} />
      ))} */}
      {/* New RS plot */}
      <RsProfileFormPlot />
    </div>
  );
}
