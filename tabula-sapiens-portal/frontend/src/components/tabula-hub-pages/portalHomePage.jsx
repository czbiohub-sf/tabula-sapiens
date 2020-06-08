import React, { Component } from "react";

import { Banner, Link as CZUILink } from "cz-ui";

import InfoBox from "~/components/basic/InfoBox.jsx";
import Heading from "~/components/basic/Heading.jsx";
import cs from "./portalHomePage.module.scss";

class HomePage extends Component {
  render() {
    return (
      <div>
        <Banner
          backgroundUrl={"../../images/tabula_hub.png"}
          mainText="TABULA HUB"
          paragraph="A catalog of whole animal single cell transcriptomes"
        />
        <div className={cs.content}>
          <Heading title="TABULA SAPIENS" />
          <InfoBox
            title="Tabula Sapiens"
            description="First-draft human cell atlas of two million cells from 25 organs of eight normal human subjects"
            buttonTitle="View the Project"
            buttonLink="/sapiens"
            image={"../../images/sapiens_logo.png"}
          />
        </div>
        <div className={cs.content}>
          <Heading title="Tabula Muris" />
          <InfoBox
            title={"Tabula Muris & Tabula Muris Senis"}
            description={
              "Comprehensive resources for the cell biology community offering a detailed molecular and cell-type specific portrait of mouse transcriptomics across the animal lifespan."
            }
            buttonTitle="View The Project"
            buttonLink="https://tabula-muris-senis.ds.czbiohub.org/"
            image={"../../images/TM_logo.png"}
          />
        </div>
        <div className={cs.content}>
          <Heading title="Tabula Microcebus" />
          <InfoBox
            title="Tabula Microcebus"
            description="The mouse lemur single cell transcriptomic atlas is the first whole animal cell atlas of a non-primate organism"
            buttonTitle="View The Project"
            buttonLink=""
            target="_blank"
            image={"../../images/lemur_logo.png"}
          />
        </div>
        <div className={cs.content}>
          <Heading title="Tabula Rhinolophidae" />
          <InfoBox
            title="Tabula Rhinolophidae"
            description="The COVID-19 Bay Area Serology Testing Consortium is dedicated to 
           developing serological assays and rapid dissemination 
           of knowledge derived from testing and validation."
            buttonTitle="View The Project"
            buttonLink=""
            target="_blank"
            image={"../../images/bat_logo.png"}
          />
        </div>
      </div>
    );
  }
}

export default HomePage;
