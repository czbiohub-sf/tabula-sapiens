import React, { Component } from "react";
// import { orderBy, filter } from "lodash/fp";

import { Heading, Search, Select } from "cz-ui";
import cs from "./CellProfiles.module.scss";

// import Search from "../basic/Search.jsx";
// import Select from "../basic/Select.jsx";
// import Table from "./Table.jsx";
// import cs from "./TableStory.module.scss";

const data = [
  {
    name: "Rabia Avila",
    status: "Trained",
    location: "Lung",
    image:
      "https://www.belr.com/wp-content/uploads/2017/06/avatar-placeholder-generic-1.jpg",
  },
  {
    name: "Arlo Morrison",
    status: "Veteran",
    location: "Pancreas",
    image:
      "https://www.belr.com/wp-content/uploads/2017/06/avatar-placeholder-generic-1.jpg",
  },
  {
    name: "Rio Watkins",
    status: "Trained",
    location: "Muscle",
    image:
      "https://www.belr.com/wp-content/uploads/2017/06/avatar-placeholder-generic-1.jpg",
  },
  {
    name: "Romario Waters",
    status: "Shadowed",
    location: "Liver",
    image:
      "https://www.belr.com/wp-content/uploads/2017/06/avatar-placeholder-generic-1.jpg",
  },
  {
    name: "Judith Hulme",
    status: "Observed",
    location: "Bladder",
    image:
      "https://www.belr.com/wp-content/uploads/2017/06/avatar-placeholder-generic-1.jpg",
  },
];

const locations = ["All Tissues"].concat(
  data
    .map((person) => person.location)
    .filter((location, i, self) => self.indexOf(location) === i)
);

const statuses = [
  "All Compartments",
  "Endothelial",
  "Epithelial",
  "Stromal",
  "Immune",
];

const CellProfiles = () => {
  const [searchedValues, setSearchedValues] = React.useState([]);
  const [selectedLocationValue, setSelectedLocationValue] = React.useState(0);
  const [selectedStatusValue, setSelectedStatusValue] = React.useState(0);

  let filteredData = data
    .filter((person) => {
      return (
        (selectedLocationValue === 0 ||
          person.location === locations[selectedLocationValue]) &&
        (selectedStatusValue === 0 ||
          person.status === statuses[selectedStatusValue])
      );
    })
    .filter((person) => {
      return (
        !searchedValues ||
        searchedValues.length === 0 ||
        searchedValues.some((searchValue) => {
          return Object.values(person).some((d) => {
            return d.toLowerCase().search(searchValue) >= 0;
          });
        })
      );
    });

  const handleSearch = (event) => {
    setSearchedValues(
      event.target.value
        .split(" ")
        .filter(() => true)
        .map((s) => s.toLowerCase())
    );
  };

  const handleLocationChange = (event) => {
    setSelectedLocationValue(event.target.value);
  };

  const handleStatusChange = (event) => {
    setSelectedStatusValue(event.target.value);
  };

  return (
    <div>
      <Heading className={cs.heading} size="xl">
        Cell Profiles
      </Heading>

      <div className={cs.selectorsPane}>
        <div className={cs.search}>
          <Search onSearch={handleSearch} />
        </div>
        <div className={cs.select}>
          <Select
            items={locations.map((location, i) => ({
              value: i,
              label: location,
            }))}
            onChange={handleLocationChange}
            value={selectedLocationValue}
          />
        </div>
        <div className={cs.select}>
          <Select
            items={statuses.map((status, i) => ({ value: i, label: status }))}
            onChange={handleStatusChange}
            value={selectedStatusValue}
          />
        </div>
      </div>

      <img src={"../../../images/cellprofiles.png"} width="1600" />
    </div>
  );
};

export default CellProfiles;
