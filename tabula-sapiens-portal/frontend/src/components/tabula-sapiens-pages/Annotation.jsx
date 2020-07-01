import React, { Component } from "react";
import axios from "axios";

import Container from "@material-ui/core/Container";
import Box from "@material-ui/core/Box";
import Typography from "@material-ui/core/Typography";

import Heading from "./Heading.jsx";

// import { BiohubUI, Button, Link as CZUILink } from "cz-ui";
import Fab from "@material-ui/core/Fab";
import AddIcon from "@material-ui/icons/Add";
import Button from "@material-ui/core/Button";
import CloudUploadIcon from "@material-ui/icons/CloudUpload";

export default class AnnotationsUpdate extends Component {
  constructor(props) {
    super(props);

    this.onFileChange = this.onFileChange.bind(this);
    this.onSubmit = this.onSubmit.bind(this);

    this.state = {
      profileImg: "",
    };
  }

  onFileChange(e) {
    this.setState({ profileImg: e.target.files[0] });
  }

  onSubmit(e) {
    e.preventDefault();
    const formData = new FormData();
    formData.append("profileImg", this.state.profileImg);
    axios
      .post("http://localhost:3000/api/user-profile", formData, {})
      .then((res) => {
        console.log(res);
      });
  }

  render() {
    return (
      <div>
        <Heading title="How to reannotate the dataset" />
        <Container maxWidth="lg">
          <Box my={4}>
            {/* <Typography
              variant="h6"
              component="h3"
              gutterBottom
              align="justify"
            >
              Add your annotations here <p />
              Layout very confusing
              <p />
              We are getting there
            </Typography> */}

            <div>
              <img src={"../../images/annotations.png"} width="1000" />
            </div>
          </Box>
        </Container>

        <Heading title="Upload annotations" />

        <Container maxWidth="lg">
          <Box my={4}>
            <div className="container">
              <div className="row">
                <form onSubmit={this.onSubmit}>
                  <div className="form-group">
                    {/* <input type="file" onChange={this.onFileChange} /> */}

                    <input
                      ref={"file-upload"}
                      type="file"
                      onChange={this.onFileChange}
                    />
                    <Button
                      onClick={(e) => {
                        this.refs["file-upload"].click();
                      }}
                    ></Button>
                  </div>

                  <div className="form-group">
                    {/* <button className="btn btn-primary" type="submit">
                Upload
              </button> */}
                    <Button
                      variant="contained"
                      color="default"
                      startIcon={<CloudUploadIcon />}
                    >
                      Upload
                    </Button>
                  </div>
                </form>
              </div>
            </div>
          </Box>
        </Container>
      </div>
    );
  }
}
