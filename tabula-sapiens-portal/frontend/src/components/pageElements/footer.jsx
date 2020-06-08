import React, { Component } from "react";

import {
  NavBar,
  Link as CZUILink,
} from "cz-ui";

class Footer extends Component {
  constructor(props) {
    super(props);
    this.state = {
      // selectedIndex: 0,
    };
  }
  render() {
    return (
      <NavBar
        title={"Chan Zuckerberg Biohub (2020)"}
      />
    );
  }
}


export default Footer;
