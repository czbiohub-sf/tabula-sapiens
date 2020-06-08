import React, { Component } from "react";
import { BrowserRouter as Router, Route, Link, Switch } from "react-router-dom";

import { BiohubUI, Banner, NavBar, Link as CZUILink } from "cz-ui";

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
                  to="/protocol"
                  onClick={() => this.setState({ selectedIndex: 0 })}
                >
                  PROTOCOLS
                </CZUILink>,
                <CZUILink
                  component={Link}
                  to="/pilot1"
                  onClick={() => this.setState({ selectedIndex: 1 })}
                >
                  TSP1
                </CZUILink>,
                <CZUILink
                  component={Link}
                  to="/pilot2"
                  onClick={() => this.setState({ selectedIndex: 2 })}
                >
                  TSP2
                </CZUILink>,
              ]}
              navSelectedLinkIndex={this.state.selectedIndex}
            />
            <Switch>
              <Route path="/protocol">
                <div>
                  <img
                    src="../../images/sapiens_logo.png"
                    alt="TSP"
                    width="100%"
                  />
                </div>
              </Route>
              <Route path="/pilot1">WIP</Route>
              <Route path="/pilot2">WIP</Route>
              <Route path="/">
                <Banner
                  backgroundUrl={"../../images/sapiens_logo.png"}
                  mainText="Tabula Sapiens"
                  paragraph="Welcome to the most exciting Tabula so far!"
                />
              </Route>
            </Switch>
          </div>
        </Router>
      </BiohubUI>
    );
  }
}

export default AppSAPIENS;
