import React from 'react';
import { Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import ReferenceSignatures from './referenceSignatures';
import MutationalSigPro from './mutationalSignatureProfile';
import { dispatchError, dispatchExploring } from '../../../services/store';

const { Header, Body } = Card;
const { Item, Link } = Nav;

export default function Explore() {
  const rootURL = window.location.pathname;

  const { displayTab, projectID } = useSelector((state) => state.exploring);

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
              { title: 'Reference Signatures', id: 'referenceSignatures' },
              { title: 'Mutational Signature Profile', id: 'mutationalSigPro' },
              // { title: 'Cosine Similarity', id: 'cosineSimilarity' },
              // {
              //   title: 'Mutational Pattern Enrichment Analysis',
              //   id: 'mutationalPattern',
              // },
              // { title: 'Profile Comparison', id: 'profileComparison' },
              // { title: 'PCA', id: 'pca' },
              // { title: 'Download', id: 'download' },
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
            display: displayTab == 'referenceSignatures' ? 'block' : 'none',
          }}
        >
          <ReferenceSignatures submitR={(fn, args) => submitR(fn, args)} />
        </Body>
        <Body
          style={{
            display: displayTab == 'mutationalSigPro' ? 'block' : 'none',
          }}
        >
          <MutationalSigPro submitR={(fn, args) => submitR(fn, args)} />
        </Body>
        {/*  <Body
          style={{
            display: displayTab == 'cosineSimilarity' ? 'block' : 'none',
          }}
        >
          <CosineSimilarity
            getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
            downloadResults={(path) => downloadResults(path)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        </Body>
        <Body
          style={{
            display: displayTab == 'mutationalPattern' ? 'block' : 'none',
          }}
        >
          <MutationalPattern
            downloadResults={(path) => downloadResults(path)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        </Body>

        <Body
          style={{
            display: displayTab == 'profileComparison' ? 'block' : 'none',
          }}
        >
          <ProfileComparison
            getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        </Body>
        <Body style={{ display: displayTab == 'pca' ? 'block' : 'none' }}>
          <PCA
            getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
            downloadResults={(path) => downloadResults(path)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        </Body>
        <Body style={{ display: displayTab == 'download' ? 'block' : 'none' }}>
          <Download />
        </Body> */}
      </Card>
    </div>
  );
}
