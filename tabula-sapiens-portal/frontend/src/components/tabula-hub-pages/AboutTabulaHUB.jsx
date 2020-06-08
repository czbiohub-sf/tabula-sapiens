import React, { Component } from "react";
import Container from "@material-ui/core/Container";
import Box from "@material-ui/core/Box";
import Typography from "@material-ui/core/Typography";

class AboutTabulaHUB extends Component {
  render() {
    return (
      <div>
        <Container maxWidth="lg">
          <Box my={4}>
            <Typography
              variant="h6"
              component="h1"
              gutterBottom
              align="justify"
            >
              The motivation behind the different Tabula projects is to create a
              set of comprehensive resources for the cell biology community
              which offers a detailed molecular and cell-type specific portrait
              of the different animals. We view such cell atlas as an essential
              companion to the genome: the genome provides a blueprint for the
              organism but does not explain how genes are used in a cell type
              specific manner or how the usage of genes changes over the
              lifetime of the organism. The cell atlas provides a deep
              characterization of phenotype and physiology which can serve as a
              reference for understanding many aspects of the cell biological
              changes that mammals undergo during their lifespan. All these are
              resources for the community and therefore free of use! The data
              and code for each project are available and we continously support
              their use. One critical factor in the Tabula projects was our
              large collaborative network of PI’s with deep expertise at
              preparation of diverse organs, enabling all organs from a subject
              to be successfully processed within a single day. We have built
              the logistics and infrastructure capable of tracking hundreds of
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
      </div>
    );
  }
}

export default AboutTabulaHUB;
