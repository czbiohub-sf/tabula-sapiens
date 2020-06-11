import React, { Component } from "react";
import { BrowserRouter as Router, Route, Link, Switch } from "react-router-dom";

import { BiohubUI, Banner, NavBar, Link as CZUILink } from "cz-ui";

import AboutSapiens from "./tabula-sapiens-pages/About.jsx";
import SapiensWorkflows from "./tabula-sapiens-pages/Workflows.jsx";
import CellProfiles from "./tabula-sapiens-pages/CellProfiles.jsx";

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
        <Router basename="/sapiens">
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
                  to="/annotation"
                  onClick={() => this.setState({ selectedIndex: 3 })}
                >
                  CELL TYPE ANNOTATION
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
                <CZUILink
                  component={Link}
                  to="/immune"
                  onClick={() => this.setState({ selectedIndex: 6 })}
                >
                  IMMUNE COMPARTMENT
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
                <img src={"../../images/annotations.png"} />
              </Route>
              <Route path="/tissues">WIP</Route>
              <Route path="/cells">
                <CellProfiles />
              </Route>
              <Route path="/immune">
                <img src={"../../images/immune.png"} width="1600" />
              </Route>
              <Route path="/home">
                <Banner
                  backgroundUrl={"../../images/sapiens_logo.png"}
                  mainText="Tabula Sapiens"
                  paragraph="Welcome to the most exciting Tabula so far!"
                />
                <img src={"../../images/dashboard.png"} />
              </Route>
            </Switch>
          </div>
        </Router>
      </BiohubUI>
    );
  }
}

export default AppSAPIENS;
