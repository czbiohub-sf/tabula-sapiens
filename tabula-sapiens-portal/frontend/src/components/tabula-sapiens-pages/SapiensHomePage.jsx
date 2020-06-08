import React, { Component } from "react";
import Container from "@material-ui/core/Container";
import Box from "@material-ui/core/Box";
import Typography from "@material-ui/core/Typography";

import { BiohubUI, Button, Link as CZUILink } from "cz-ui";

import AppSAPIENS from "../appSAPIENS.jsx";

class SapiensHomePage extends Component {
  render() {
    return (
      <div>
        <AppSAPIENS />
      </div>
    );
  }
}

export default SapiensHomePage;
