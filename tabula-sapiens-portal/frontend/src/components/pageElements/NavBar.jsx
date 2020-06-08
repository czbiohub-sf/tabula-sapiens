import React, { Component } from "react";
import { BrowserRouter as Router, Link } from "react-router-dom";

import { Button, NavBar, Link as CZUILink } from "cz-ui";

class NavBarPortal extends Component {
  constructor(props) {
    super(props);
    this.state = {
      selectedIndex: 0,
    };
  }
  render() {
    return (
      <NavBar
        accent
        title={"Tabula HUB"}
        titleOptions={[
          <CZUILink
            component={Button}
            href="https://tabula-muris.ds.czbiohub.org/"
            target="_blank"
          >
            Tabula Muris
          </CZUILink>,
          <CZUILink component={Button} 
          href="https://tabula-muris-senis.ds.czbiohub.org/"
          target="_blank">
            Tabula Muris Senis
          </CZUILink>,
          <CZUILink
            component={Button}
            href="http://tabula-sapiens.ds.czbiohub.org/"
            target="_blank"
          >
            Tabula Sapiens
          </CZUILink>,
        ]}
        username={" "}
        navLinks={[
          <CZUILink
            component={Link}
            to="/"
            onClick={() => this.setState({ selectedIndex: 0 })}
          >
            TABULA HUB
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
            to="/publications"
            onClick={() => this.setState({ selectedIndex: 2 })}
          >
            PUBLICATIONS
          </CZUILink>,
        ]}
        navSelectedLinkIndex={this.state.selectedIndex}
      />
    );
  }
}

export default NavBarPortal;
