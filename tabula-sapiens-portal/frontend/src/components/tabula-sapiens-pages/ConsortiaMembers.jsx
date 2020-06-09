import React from "react";
import { Grid } from "@material-ui/core";

import cs from "./ConsortiaMembers.scss";

const ConsortiaMemberBoxes = function () {
  return (
    <Grid className={cs.container}>
      <Grid item xs={3}>
        <img className={cs.biohub} src={"../../../images/biohublogo.png"} />
      </Grid>
      <Grid item xs={3}>
        <img className={cs.ucsf} src={"../../../images/ucsf_logo.jpg"} />
      </Grid>
      <Grid item xs={3}>
        <img
          className={cs.berkeley}
          src={"../../../images/UCBerkeley_logo.jpg"}
        />
      </Grid>
      <Grid item xs={3}>
        <img
          className={cs.stanford}
          src={"../../../images/stanford_logo.png"}
        />
      </Grid>
    </Grid>
  );
};

export default ConsortiaMemberBoxes;
