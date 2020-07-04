import React from "react";
import ReactDOM from "react-dom";
import AppSAPIENS from "./components/appSAPIENS.jsx";
import { BrowserRouter } from "react-router-dom";

ReactDOM.render(
  <BrowserRouter>
    <AppSAPIENS />
  </BrowserRouter>,

  document.getElementById("root")
);
