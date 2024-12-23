# Methylation-networks

**This repository contains the code and pre-processed data used in the study of tissue- and time-specific DNA methylation dynamics in Sus scrofa via joint fused ridge network modelling.**

---

## Table of Contents
1. [Usage](#usage)
2. [Contact](#contact)
3. [License](#license)

---

## Usage
To clone the repo use `git clone https://github.com/wachalak/Methylation-networks.git`

The **Code** directory contains two files: pre_processing.sh and analyses.R

- pre_processing.sh gives the workflow to obtain bed files needed for the analyses with R. The initial data used in this study were collected by the GENE-SWitCH (https://www.gene-switch.eu/), and are available through the FAANG data portal at https://data.faang.org/projects/GENE-SWitCH. The accession code is: (i) PRJEB41822 for
RRBS, (ii) PRJEB41970 for RNA-Seq, and (iii) PREJEB70458 for ChIP-Seq
- analyses.R gives the complete workflow for joint fused ridge network modelling and analysis

The **Data-post-processing** directory contains all files needed to run analyses.R. The bed files have been compressed. To decompress all bed files in the directory you can run either:

- `gunzip *.bed.gz`, or
- `for file in *.bed.gz; do bgzip -d "$file"; done`

---

## Contact

For questions or feedback, please contact:

- **Author:** Karolina Wachala  
- **Email:** karolina.wachala@unibe.ch  

---

## License

This repository is licensed under the **MIT License**. See the [LICENSE](./LICENSE) file for more details.

---

