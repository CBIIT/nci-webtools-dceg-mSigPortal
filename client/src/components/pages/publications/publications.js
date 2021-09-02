import React, { useMemo } from 'react';
import { useSelector, useDispatch } from 'react-redux';
import Table from '../../../components/controls/table/table';
import { actions } from '../../../services/store/publications';
import './publications.scss';

export default function Publications() {
  // data is retrieved in components/app.js
  const dispatch = useDispatch();
  const { orA, orB, rp, cm } = useSelector((state) => state.publications);
  const {
    globalFilter: oraSearch,
    hidden: oraHidden,
    pagination: oraPagination,
    ...oraTable
  } = orA;
  const {
    globalFilter: orbSearch,
    hidden: orbHidden,
    pagination: orbPagination,
    ...orbTable
  } = orB;
  const {
    globalFilter: rpSearch,
    hidden: rpHidden,
    pagination: rpPagination,
    ...rpTable
  } = rp;
  const {
    globalFilter: cmSearch,
    hidden: cmHidden,
    pagination: cmPagination,
    ...cmTable
  } = cm;
  const oraMemo = useMemo(() => oraTable, [oraTable]) || {};
  const orbMemo = useMemo(() => orbTable, [orbTable]) || {};
  const rpMemo = useMemo(() => rpTable, [rpTable]) || {};
  const cmMemo = useMemo(() => cmTable, [cmTable]) || {};

  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <div className="mb-4">
          <p>
            An overview of published papers, tools, websites and databases
            related to mutational signatures analysis.
          </p>
        </div>

        {oraMemo.data && (
          <div className="mb-5">
            <Table
              title="Original Research Papers Including Specific Mutational Signatures in mSigPortal"
              columns={oraMemo.columns}
              data={oraMemo.data}
              hidden={oraHidden}
              globalFilter={oraSearch}
              pagination={oraPagination}
              mergeState={(state) =>
                dispatch(actions.mergeState({ orA: state }))
              }
            />
          </div>
        )}
        {orbMemo.data && (
          <div className="mb-5">
            <Table
              title="Original Research Papers involved Mutational Signature Analyses"
              columns={orbMemo.columns}
              data={orbMemo.data}
              hidden={orbHidden}
              globalFilter={orbSearch}
              pagination={orbPagination}
              mergeState={(state) =>
                dispatch(actions.mergeState({ orB: state }))
              }
            />
          </div>
        )}
        {rpMemo.data && (
          <div className="mb-5">
            <Table
              title="Review Papers Focued on Mutational Signatures"
              columns={rpMemo.columns}
              data={rpMemo.data}
              hidden={rpHidden}
              globalFilter={rpSearch}
              pagination={rpPagination}
              mergeState={(state) =>
                dispatch(actions.mergeState({ rp: state }))
              }
            />
          </div>
        )}
        {cmMemo.data && (
          <div className="mb-5">
            <Table
              title="Computational Methods, Databases or Websites for Mutational Signature Analyses"
              columns={cmMemo.columns}
              data={cmMemo.data}
              hidden={cmHidden}
              globalFilter={cmSearch}
              pagination={cmPagination}
              mergeState={(state) =>
                dispatch(actions.mergeState({ cm: state }))
              }
            />
          </div>
        )}

        <div className="mb-4">
          <h3>Citations</h3>
        </div>

        <p>Last update: 20 JAN 2021.</p>
      </div>
    </div>
  );
}
