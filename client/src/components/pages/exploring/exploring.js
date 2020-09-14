import React, { useEffect } from 'react';
import { Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import SignatureExploring from './signatureExploring';
import {
  dispatchError,
  dispatchExploring,
  dispatchExpMutationalProfiles,
} from '../../../services/store';

const { Header, Body } = Card;
const { Item, Link } = Nav;

export default function Explore() {
  const rootURL = window.location.pathname;

  const { displayTab, projectID, referenceSignatureData } = useSelector(
    (state) => state.exploring
  );

  useEffect(async () => {
    // if (!referenceSignatureData.size) getReferenceSignatures();
  }, []);

  async function getReferenceSignatures() {
    const data = await (
      await fetch(`${rootURL}getReferenceSignatures`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
      })
    ).json();

    const signatureSourceOptions = [...new Set(data.map((row) => row.Source))];
    const signatureSource = signatureSourceOptions[0];
    let filteredData = data.filter((row) => row.Source == signatureSource);
    const profileNameOptions = [
      ...new Set(filteredData.map((row) => row.Profile)),
    ];
    const profileName = profileNameOptions[0];
    filteredData.filter(
      (row) => row.Source == signatureSource && row.Profile == profileName
    );
    const refSignatureSetOptions = [
      ...new Set(filteredData.map((row) => row.Signature_set_name)),
    ];
    const refSignatureSet = refSignatureSetOptions[0];
    filteredData.filter(
      (row) =>
        row.Source == signatureSource &&
        row.Profile == profileName &&
        row.Signature_set_name == refSignatureSet
    );
    const strategyOptions = [
      ...new Set(filteredData.map((row) => row.Dataset)),
    ];
    const strategy = strategyOptions[0];
    filteredData.filter(
      (row) =>
        row.Source == signatureSource &&
        row.Profile == profileName &&
        row.Signature_set_name == refSignatureSet &&
        row.Dataset == strategy
    );
    const signatureNameOptions = [
      ...new Set(filteredData.map((row) => row.Signature_name)),
    ];

    dispatchExploring({
      referenceSignatureData: data,
    });
    dispatchExpMutationalProfiles({
      signatureSource: signatureSource,
      signatureSourceOptions: signatureSourceOptions,
      profileName: profileName,
      profileNameOptions: profileNameOptions,
      refSignatureSet: refSignatureSet,
      refSignatureSetOptions: refSignatureSetOptions,
      strategy: strategy,
      strategyOptions: strategyOptions,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    });
  }

  function submitR(fn, args) {
    return fetch(`${rootURL}exploringR`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ fn: fn, args: args, projectID: projectID }),
    });
  }

  return (
    <div className="position-relative m-3">
      <Card>
        <Header>
          <Nav variant="pills" defaultActiveKey="#mutationalProfiles">
            {[
              { title: 'Signature Exploring', id: 'signatureExploring' },
              { title: 'Exposure Exploring', id: 'exposureExploring' },
            ].map(({ title, id }) => {
              return (
                <Item key={id}>
                  <Link
                    active={displayTab == id}
                    onClick={() => dispatchExploring({ displayTab: id })}
                  >
                    {title}
                  </Link>
                </Item>
              );
            })}
          </Nav>
        </Header>
        <Body
          style={{
            display: displayTab == 'signatureExploring' ? 'block' : 'none',
          }}
        >
          <SignatureExploring />
        </Body>
        {/* <Body
          style={{
            display: displayTab == 'exposureExploring' ? 'block' : 'none',
          }}
        >
          <MutationalProfiles />
        </Body> */}
      </Card>
    </div>
  );
}
