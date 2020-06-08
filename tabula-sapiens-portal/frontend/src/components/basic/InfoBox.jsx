import React from "react";

import cs from "./InfoBox.scss";

const InfoBox = function (props) {
  return (
    <div className={cs.container}>
      <img className={cs.image} src={props.image}></img>
      <div className={cs.text}>
        <h1 className={cs.title}>{props.title}</h1>
        <p>{props.description}</p>
        <a href={props.buttonLink} className={cs.buttonLink}>
          <div className={cs.button}>{props.buttonTitle}</div>
        </a>
      </div>
    </div>
  );
};

export default InfoBox;
