import { atom, selector } from 'recoil';
import axios from 'axios';

export const defaultTreeLeafData = { links: [], nodes: [] };

export const getGraphData = selector({
  key: 'treeLeaf.plotData',
  get: async ({ get }) => {
    try {
      const { data } = await axios.post('api/visualizationWrapper', {
        fn: 'getTreeLeaf',
      });

      return data.output;
    } catch (error) {
      return null;
    }
  },
});
