import { memo, useEffect, useState } from 'react';

function TabPanel({ active, children }) {
  const [visited, setVisited] = useState(active);
  if (active && !visited) setVisited(true);

  // plots can't measure while hidden, so let them refit when shown
  useEffect(() => {
    if (active) window.dispatchEvent(new Event('resize'));
  }, [active]);

  if (!visited) return null;
  return <div className={active ? 'd-block' : 'd-none'}>{children}</div>;
}

// hidden tabs skip re-renders and catch up with the latest props when shown
export default memo(TabPanel, (prev, next) => !prev.active && !next.active);
