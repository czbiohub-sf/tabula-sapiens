import React , {useState, useEffect} from 'react';
import { useForm } from 'react-hook-form'
import './App.css';
import {
    BiohubUI,
    AppBar,
    Select,
    Input,
} from "cz-ui";
import cs from "./App.module.scss";
import { Button } from '@material-ui/core';

const ann_methods = [
    "All Methods",
    "OnClass",
    "SCANVI",
    "SVM",
    "SingleCellNet",
];
const use_gpu = [
    "Yes",
    "No",
];

const tissues = [
    'All Tissues',
    "Bladder",
    "Blood",
    "Bone Marrow",
    "Kidney",
    "Large Intestines",
    "Lung",
    "Lymph Node",
    "Muscle",
    "Pancreas",
    "Skin",
    "Small Intestine",
    "Spleen",
    "Thymus",
    "Trachea",
    'Vasculature'
]

const use_10x_only = [
    'True',
    'False'
]

const batch_correction_conditions = [
    'No Batch Correction',
    'Donor',
    'Technology'

]
function App() {
    const [selectedLocationValue, setSelectedLocationValue] = React.useState(0);
    const handleLocationChange = (event) => {
        setSelectedLocationValue(event.target.value);
    };
    const { register, handleSubmit } = useForm()

    const onSubmit = async (data) => {
        const formData = new FormData()
        formData.append("adata", data.adata[0])

        const res = await fetch("http://localhost:4000/adata", {
            method: "POST",
            body: formData
        }).then(res => res.json())
        alert(JSON.stringify(res))
    }
    const [selectedStatusValue, setSelectedStatusValue] = React.useState(0);
    const handleStatusChange = (event) => {
        setSelectedStatusValue(event.target.value);
    };

    const [currentTime, setCurrentTime] = useState(0);

    useEffect(() => {
        fetch('/time').then(res => res.json()).then(data => {
          setCurrentTime(data.time);
        });
    }, []);

    return (
	    <BiohubUI>
            <div role="main">
                <AppBar title={"czi.team"} rightHeader={"Automated Annotation with Tabula Sapiens"} />
                <form onSubmit={handleSubmit(onSubmit)}>
                    <div className={cs.input}>
                        Input Anndata:
                    <input ref={register} type="file" name="adata" />
                        <Button ref={register} type={"hollow-dark"}>Submit AnnData</Button>
                    </div>
                </form>
                <div className={cs.selectorsPane}>
                    <div className={cs.txt}>
                        Annotation Method:
                    </div>
                    <div className={cs.select}>
                        <Select
                            items={ann_methods.map((ann_method, i) => ({ value: i, label: ann_method }))}
                            onChange={handleLocationChange}
                            value={selectedLocationValue}

                        />
                    </div>
                </div>
        	<p>The current time is {currentTime}.</p>
                <div className={cs.selectorsPane}>

                    <div className={cs.txt}>
                        Tissue:
                    </div>
                    <div className={cs.select}>
                        <Select
                            items={tissues.map((status, i) => ({ value: i, label: status }))}
                            onChange={handleLocationChange}
                            value={selectedLocationValue}
                        />
                    </div>
                    <div className={cs.txt}>
                        Batch Correction Conditions:
                    </div>
                    <div className={cs.select}>
                        <Select
                            items={batch_correction_conditions.map((status, i) => ({ value: i, label: status }))}
                            onChange={handleLocationChange}
                            value={selectedLocationValue}
                        />
                    </div>
                </div>
                <div className={cs.selectorsPane}>
                    <div className={cs.txt}>
                        Use 10X Data Only:
                    </div>
                    <div className={cs.select}>
                        <Select
                            items={use_10x_only.map((status, i) => ({ value: i, label: status }))}
                            onChange={handleLocationChange}
                            value={selectedLocationValue}
                        />
                    </div>
                    <div className={cs.txt}>
                        Use GPU:
                    </div>
                    <div className={cs.select}>
                        <Select
                            items={use_gpu.map((status, i) => ({ value: i, label: status }))}
                            onChange={handleLocationChange}
                            value={selectedLocationValue}
                        />
                    </div>

                </div>
                <div className={cs.input}>
                    Output Folder: <Input>Testing</Input>
                </div>
                <div className={cs.input}>
                    <Button type={'hollow-dark'}>Run Annotation Pipeline</Button>
                </div>
            </div>
        </BiohubUI>
    );
}

export default App;
