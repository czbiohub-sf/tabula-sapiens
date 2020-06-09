import React, { Component } from "react";
import Container from "@material-ui/core/Container";
import Box from "@material-ui/core/Box";
import Typography from "@material-ui/core/Typography";

import Heading from "./Heading.jsx";
import ConsortiaMemberBoxes from "./ConsortiaMembers.jsx";

class AboutSapiens extends Component {
  render() {
    return (
      <div>
        <Heading title="About" />
        <Container maxWidth="lg">
          <Box my={4}>
            <Typography
              variant="h6"
              component="h3"
              gutterBottom
              align="justify"
            >
              Tabula Sapiens will be a benchmark, first-draft human cell atlas
              of two million cells from 25 organs of eight normal human
              subjects. <p />
              The goal of the project is to deliver annotated and analyzed
              high-quality transcriptomic data that will be integrated in the
              data coordination platform of the Human Cell Atlas. We will build
              directly on our unique skills, experience, and data infrastructure
              from Tabula Muris, Tabula Muris Senis and Tabula Microcebus to
              create a high-quality human reference dataset and portal at a
              10-fold larger scale from these prior efforts.
              <p />
              One critical factor in the Tabula projects was our large
              collaborative network of PI’s with deep expertise at preparation
              of diverse organs, enabling all organs from a subject to be
              successfully processed within a single day. We have built the
              logistics and infrastructure capable of tracking hundreds of
              samples and thousands of 384-well plates from tissue through
              sample prep, library construction and on to sequencing and
              ultimately computational and expert cell annotation with tight
              quality control. We will supplement our network with human tissue
              expertise and use our experience to balance and assign cell types
              from each tissue compartment and optimally mix high-quality
              plate-seq data and high-volume droplet-based data to provide a
              broad and deep benchmark atlas.
            </Typography>
          </Box>
        </Container>

        <Heading title="Consortium Members" />
        <ConsortiaMemberBoxes />
      </div>
    );
  }
}

export default AboutSapiens;
