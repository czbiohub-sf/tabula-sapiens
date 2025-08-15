[![DOI](https://zenodo.org/badge/200055743.svg)](https://zenodo.org/badge/latestdoi/200055743)


# tabula-sapiens

Welcome to Tabula Sapiens!

<img src="https://github.com/czbiohub/tabula-sapiens/blob/master/sapiens_logo.png" width="50%" height="50%">

# Where is the data
## Raw data
Since April 2021, [Tabula Sapiens data]("https://aws.amazon.com/marketplace/pp/prodview-3knq6mp7oquu6") have been made available to all users free of charge. This product is part of the AWS Open Data Sponsorship Program and contains data sets that are publicly available for anyone to access and use. AWS has made the data freely available so that anyone can download the resource to perform analysis and advance medical discovery without needing to worry about the cost of storing Tabula Sapiens data or the time required to download it. The data can be browsed but before you can download we require users to [complete a data sharing agreement](https://apps.docusign.com/webforms/us/929255c461f4bb15efc65d776f103035?r=1).

The s3 public bucket follows this example structure:

```
aws s3 ls czb-tabula-sapiens
├── Donor1/
├── Donor10/
├── Donor11/
├── Donor12/
├── Donor13/
├── Donor14/
├── Donor15/
├── Donor2/
│   ├── alignment-gencode/
│   │   ├── 10X/
│   │   ├── smartseq2/
│   ├── fastqs/
│   │   ├── 10X/
│   │   ├── smartseq2/
│   ├── gene-count-tables/
│   ├── immune-repertoire-analysis/
│   │   ├── bracer/
│   │   ├── tracer/
├── Donor3/
├── Donor4/
├── Donor5/
├── Donor6/
├── Donor7/
├── Donor8/
├── Donor9/
├── reference/
│   ├── cellranger/
│   │   ├── homo.gencode.v30.annotation.ERCC92.tgz
│   ├── STAR/
│   │   ├── homo.gencode.v30.annotation.ERCC92.tgz

```

Each `DonorN` folder contains all the raw data for `TSPN`. Each `fastqs` folder is split between `10x` and `smartseq` with respective files inside, an identical strucutre followed by the `alignment-gencode` folder. The `gene-count-tables` correspond the smartseq2 only. For the donors for which we generated smartseq2 we also provide the raw outputs for `bracer` and `tracer` inside the respective `immune-repertoire-analysis` folder. The reference files to use to reprocess the entire dataset from fastqs onwards are provide in `reference`.

For instructions on how to create an AWS account (free of charge) please refer to [AWS documentation](https://aws.amazon.com/premiumsupport/knowledge-center/create-and-activate-aws-account/)


## Processed data
Our ready-to-use data is available from figshare: https://figshare.com/projects/Tabula_Sapiens/100973

If you would like to use cellxgene to explore the data on your local machine consider using [exploratory-cellxgene](https://github.com/czbiohub/cellxgene)



# How to use cellxgene
Checkout detailed instructions [here](run-cellxgene.md)

# Data Portal
To interact with the data checkout our portal: http://tabula-sapiens-portal.ds.czbiohub.org/
