import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import Card from 'react-bootstrap/Card';
import { CardDeck, Button } from 'react-bootstrap';
import './home.scss';



export default function Home({ links }) {

  const [cats,setCats] = useState([])

  const cat = async () => {
    const response = await fetch('http://aws.random.cat/meow')
    const obj = await response.json()
    return obj.file
  }

  useEffect(() => {
    for(let i in links.slice(1,5)){
      cat().then((url) => setCats((cats) => [...cats,url]))
    }
  },[])

  const Rows = (props) =>{
    
    var start = Number(props.start)
    var end = Number(props.end)
    var index = start-2
    console.log(start)
    
    return(
      <CardDeck>
          {links
            .slice(start, end)
            .map(
              (
                { exact, route, action, title, cardTitle, cardText, image }
              ) => {
                index = index+1
                return (
                <>
                  
                  <Card
                    key={title}
                    className="mb-5 align-self-center"
                    style={{
                      width: '18rem',
                      justifyContent: 'center',
                      alignItems: 'center',
                      border: '1px solid #DADBE6',
                      // minHeight: '280px'
                      // borderRadius: '10px'
                    }}
                  >
                    <Link
                      className="stretched-link"
                      exact={exact}
                      key={index}
                      to={route}
                    >
                    
                    <span className="ssr-only" display>{title + ' link'}</span> {/*Causes style issues if image is too small*/}
                      <div
                        className="bg-primary rounded-circle"
                        style={{ marginTop: '-40px', padding: '10px' }}
                      >
                       <img alt="icon" src={cats[index]} height="100" width="100" /> {/*Doesn't fit in circle properly*/}
 
                      </div>
                    </Link> 
                    <Card.Body>
                      <Card.Title
                        style={{ color: '#545871', wordSpacing: '100vw' }}
                      >
                        <h2 style={{ fontSize: '1.75rem' }}>
                          <b>{cardTitle}</b>
                        </h2>
                      </Card.Title>
                      <Card.Text className="text-secondary">
                        <small>{cardText}</small>
                      </Card.Text>
                    </Card.Body>
                    <Card.Footer
                      className="bg-white border-top-0"
                      style={{ width: '100%' }}
                    >
                      <Button
                        className="my-2 border border-0 font-weight-bold"
                        style={{
                          backgroundColor: '#2CC799',
                          // borderRadius: '10px',
                          width: '90%',
                        }}
                      >
                        <Link
                          className="stretched-link text-dark"
                          style={{ textDecoration: 'none' }}
                          exact={exact}
                          key={index}
                          to={route}
                        >
                          {action}
                        </Link>
                      </Button>
                    </Card.Footer>
                  </Card>
                  <div className="d-lg-none w-100"></div>
                </>
              )
              
              })}
                      
        </CardDeck>
    )
  }
  
  return (
    <>
      <div className="banner-container text-center d-none d-md-block">
        {/* <img
          src="assets/images/plco-banner.jpg"
          alt="PLCO banner"
          style={{ width: '100%' }}
        ></img> */}
        <div className="banner-overlay-text row justify-content-center text-center text-dark w-75">
          <div className="col-12">
            <h1 className="text-dark">
              <b>mSigPortal</b>
            </h1>
          </div>
          <div
            className="col-6 w-50 my-3 align-self-center"
            style={{ borderTop: '3px solid black' }}
          ></div>
          <div
            className="col-12 text-center mt-2 font-weight-bold"
            style={{ width: '100%', fontSize: '18pt' }}
          >
            App Description Line 1
            <br />
            App Description Line 2
          </div>
          <div
            className="col-12 text-center mt-5"
            style={{ width: '100%', fontSize: '14pt' }}
          >
            {/* <Button
              className="mr-5 px-4"
              style={{ backgroundColor: '#F2711D', border: 'none' }}>
              Link
            </Button>
            <Button
              className="px-4"
              style={{ backgroundColor: '#01BDD4', border: 'none' }}>
              Link
            </Button> */}
          </div>
        </div>
      </div>

      {/* mobile */}
      <div className="text-center mt-2 d-md-none">
        <h1 className="text-dark">
          <b>mSigPortal</b>
        </h1>
        <hr className="w-75"></hr>
        <div className="px-3 text-center">App Description</div>
      </div>

      <div
        className="container align-middle text-center"
        style={{ marginTop: '70px' }}
      >
      {/*Cards go here, could change to not hardcode numbers*/}
      {<Rows start={1} end={3}></Rows>}
      {<Rows start={3} end={5}></Rows>}
        
      </div>
      <div className="bg-white text-center">
        <div
          className="bg-secondary text-dark text-center"
          style={{
            height: '50px',
            clipPath: 'polygon(50% 100%, 0 0, 100% 0)',
          }}
        ></div>
        <div className="py-5">
          <h3 style={{ color: '#545871' }}>
            <b>OUR FOCUS</b>
          </h3>
          <h4 className="container mt-3 text-dark" style={{ fontSize: '16pt' }}>
            Mission Text Placeholder
          </h4>
        </div>
      </div>
      <div className="bg-egg py-4">
        <div className="container my-3 text-dark">
          Credits: TBD
          <br />
          Citation: TBD
          <br />
          mSigPortal is available under the MIT license, an Open Source
          Initiative approved license.
        </div>
      </div>
    </>
  );
}
