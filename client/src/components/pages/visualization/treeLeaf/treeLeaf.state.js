import { atom, selector } from 'recoil';
import axios from 'axios';

export const defaultTreeLeafData = { links: [], nodes: [] };

export const treeLeafData = selector({
  key: 'treeLeaf.plotState',
  get: async ({ get }) => {
    try {
      const { data } = await axios.post('api/visualizationWrapper', {
        fn: 'getTreeAndLeaf',
      });

      return data.output;
    } catch (error) {
      return defaultTreeLeafData;
    }
  },
});
