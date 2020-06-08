import React, { Component } from "react";
import ReactDOM from "react-dom";
import {
  BrowserRouter as Router,
  Route,
  Link,
  NavLink,
  Switch,
} from "react-router-dom";

import {
  BiohubUI,
  PageTitle,
  AppBar,
  Button,
  NavBar,
  Link as CZUILink,
} from "cz-ui";

import Publications from "./tabula-hub-pages/publicationsPage.jsx";
import AboutTabulaHUB from "./tabula-hub-pages/AboutTabulaHUB.jsx";
import HomePage from "./tabula-hub-pages/portalHomePage.jsx";
import SapiensHomePage from "./tabula-sapiens-pages/SapiensHomePage.jsx";

import Footer from "./pageElements/footer.jsx";
import NavBarPortal from "./pageElements/NavBar";

class App extends Component {
  render() {
    return (
      <BiohubUI>
        <Router basename="/">
          <div>
            <Switch>
              <Route path="/about">
                <NavBarPortal />
                <AboutTabulaHUB />
              </Route>
              <Route path="/publications">
                <NavBarPortal />
                <Publications />
              </Route>
              <Route path="/sapiens">
                <SapiensHomePage />
              </Route>
              <Route path="/">
                <NavBarPortal />
                <HomePage />
              </Route>
            </Switch>
            <Footer />
          </div>
        </Router>
      </BiohubUI>
    );
  }
}

export default App;
