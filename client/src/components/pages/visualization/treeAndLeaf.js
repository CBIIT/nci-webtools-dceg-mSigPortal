import axios from 'axios';
import { forceSimulation } from 'd3-force';
import { useEffect, useState } from 'react';
export default function TreeAndLeaf() {
  const [output, setOutput] = useState('');

  async function getData() {
    const { data } = await axios.post('api/visualizationWrapper', {
      fn: 'treeAndLeaf',
    });
    setOutput(JSON.stringify(data));
    console.log(data);
  }

  useEffect(() => {
    getData();
  }, []);

  return (
    <pre className="container-md">
      <code>{output}</code>
    </pre>
  );
}
