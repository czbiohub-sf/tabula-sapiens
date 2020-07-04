import React, { Component } from "react";
import { BrowserRouter as Router, Route, Link, Switch } from "react-router-dom";

import { BiohubUI, Banner, NavBar, Link as CZUILink, Footer } from "cz-ui";

import AboutSapiens from "./tabula-sapiens-pages/About.jsx";
import SapiensWorkflows from "./tabula-sapiens-pages/Workflows.jsx";
import CellProfiles from "./tabula-sapiens-pages/CellProfiles.jsx";
import AnnotationsUpdate from "./tabula-sapiens-pages/Annotation.jsx";
import MySankeyComponent from "./tabula-sapiens-pages/CellTypes.jsx";

class AppSAPIENS extends Component {
  constructor(props) {
    super(props);
    this.state = {
      selectedIndex: 0,
    };
  }
  render() {
    return (
      <BiohubUI>
        <Router basename="/">
          <div>
            <NavBar
              accent
              title={"Tabula Sapiens"}
              navLinks={[
                <CZUILink
                  component={Link}
                  to="/home"
                  onClick={() => this.setState({ selectedIndex: 0 })}
                >
                  HOME
                </CZUILink>,
                <CZUILink
                  component={Link}
                  to="/about"
                  onClick={() => this.setState({ selectedIndex: 1 })}
                >
                  ABOUT
                </CZUILink>,
                <CZUILink
                  component={Link}
                  to="/workflows"
                  onClick={() => this.setState({ selectedIndex: 2 })}
                >
                  WORKFLOWS
                </CZUILink>,
                <CZUILink
                  component={Link}
                  to="/celltypes"
                  onClick={() => this.setState({ selectedIndex: 3 })}
                >
                  CELL TYPES
                </CZUILink>,
                <CZUILink
                  component={Link}
                  to="/annotation"
                  onClick={() => this.setState({ selectedIndex: 3 })}
                >
                  CONTRIBUTE
                </CZUILink>,
                <CZUILink
                  component={Link}
                  to="/tissues"
                  onClick={() => this.setState({ selectedIndex: 4 })}
                >
                  TISSUE PROFILES
                </CZUILink>,
                <CZUILink
                  component={Link}
                  to="/cells"
                  onClick={() => this.setState({ selectedIndex: 5 })}
                >
                  CELL PROFILES
                </CZUILink>,
                // <CZUILink
                //   component={Link}
                //   to="/immune"
                //   onClick={() => this.setState({ selectedIndex: 6 })}
                // >
                //   IMMUNE COMPARTMENT
                // </CZUILink>,
                <CZUILink
                  component={Link}
                  to="/splicing"
                  onClick={() => this.setState({ selectedIndex: 7 })}
                >
                  SPLICING
                </CZUILink>,
              ]}
              navSelectedLinkIndex={this.state.selectedIndex}
            />
            <Switch>
              <Route path="/about">
                <AboutSapiens />
              </Route>
              <Route path="/workflows">
                <SapiensWorkflows />
              </Route>
              <Route path="/annotation">
                <AnnotationsUpdate />
                {/* <img src={"../../images/annotations.png"} width="1600" /> */}
              </Route>
              <Route path="/tissues">
                <img src={"../../images/tissueprofile.png"} width="1600" />
              </Route>
              <Route path="/cells">
                <CellProfiles />
              </Route>
              {/* <Route path="/immune">
                <img src={"../../images/immune.png"} width="1600" />
              </Route> */}
              <Route path="/splicing">WIP</Route>
              <Route path="/celltypes">
                <MySankeyComponent />
              </Route>
              <Route path="/">
                <Banner
                  backgroundUrl={"../../images/sapiens_logo.png"}
                  mainText="Tabula Sapiens"
                  paragraph="Welcome to the most exciting Tabula so far!"
                />
                <img src={"../../images/dashboard.png"} width="1600" />
              </Route>
            </Switch>
            <Footer navLinks={[<CZUILink href="/">HOME</CZUILink>]} />
          </div>
        </Router>
      </BiohubUI>
    );
  }
}

export default AppSAPIENS;
