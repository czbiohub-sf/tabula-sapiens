import React, { Component } from "react";
import Grid from "@material-ui/core/Grid";

import MySankeyComponent from "../sankey-components/SankeyAll.jsx";
import MySankeyEpithelialComponent from "../sankey-components/SankeyEpithelial.jsx";
import MySankeyEndothelialComponent from "../sankey-components/SankeyEndothelial.jsx";
import MySankeyStromalComponent from "../sankey-components/SankeyStromal.jsx";
import MySankeyImmuneComponent from "../sankey-components/SankeyImmune.jsx";

class CellTypesComponent extends Component {
  render() {
    return (
      <div id="all-sankeys">
        <Grid container spacing={3} justify="center" alignItems="center">
          <Grid item xs={12} alignItems="center" justify="center">
            <MySankeyComponent />
          </Grid>
          <Grid item xs={6}>
            <MySankeyEpithelialComponent />
          </Grid>
          <Grid item xs={6}>
            <MySankeyImmuneComponent />
          </Grid>
          <Grid item xs={6}>
            <MySankeyEndothelialComponent />
          </Grid>
          <Grid item xs={6}>
            <MySankeyStromalComponent />
          </Grid>
        </Grid>
      </div>
    );
  }
}

export default CellTypesComponent;
