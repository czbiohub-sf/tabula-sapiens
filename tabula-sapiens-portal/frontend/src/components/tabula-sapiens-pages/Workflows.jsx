import React, { Component } from "react";
import Container from "@material-ui/core/Container";
import Box from "@material-ui/core/Box";
import Typography from "@material-ui/core/Typography";

import Heading from "./Heading.jsx";

class SapiensWorkflows extends Component {
  render() {
    return (
      <div>
        <Heading title="GitHub" />

        <img src={"../../../images/github.png"} />

        <Heading title="protocols.io" />
        <img src={"../../../images/protocols.io.png"} />
      </div>
    );
  }
}

export default SapiensWorkflows;
