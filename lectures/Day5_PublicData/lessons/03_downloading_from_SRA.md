---
title: "Accessing public genomic data: SRA"
author: "Shahin Shahsavari"
date: "October, 2020"
---

# Downloading data from SRA

The Sequence Read Archive (SRA) is an archive for high throughput sequencing data, publically accessible, for the purpose of enhancing reproducibility in the scientific community.

There are four hierarchical levels of SRA entities and their accessions:  

1. **STUDY** with accessions in the form of SRP, ERP, or DRP  
2. **SAMPLE** with accessions in the form of SRS, ERS, or DRS  
3. **EXPERIMENT** with accessions in the form of SRX, ERX, or DRX  
4. **RUN** with accessions in the form of SRR, ERR, or DRR

The minimum publishable unit in the SRA, is an EXPERIMENT (SRX)

<img src="../img/sra_structure_infograph.png" width="600">

But most commonly, we find data we are interested in starting from a publication (or study) in GEO. 

This time let's use the other GEO dataset from the paper "GSE51443", this is the one for the "Mov10 iCLIP-SEQ" dataset. [Click here to access the GEO summary page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51443) for this second (smaller) dataset.

<img src="../img/mov10_geo.png" width="600">

Towards the bottom of the page you will find a link for **"SRA"** under the heading **"Relations"**.

<img src="../img/sra_relations.png" width="600">

Clicking on this link takes you to a page that lists all the biological samples for the study - each with a link to their specific runs and files. If we were only interested in one sample, we could follow the relevant link and find its runs. But generally we want the files for all samples and their replicates, and to find this in one comprehensive list, we use the **Run Selector**. Navigate to the bottom of the page and click **"send to"** and select **"Run Selector"**, and then press **"go"**.

<img src="../img/send_to_run_selector.png" width="600">

## Run selector
You'll notice that the Run Selector has aggregated all the information for the study samples, including a table of metadata at the top, giving information on: 

- LibraryLayout (whether the reads were sequenced using single or paired end sequencing)
- Platform (which sequencing technology was used)
- other useful information that should be noted for downstream analysis

<img src="../img/run_table.png" width="700">

Below this there is also a summary line detailing the total number of runs in the study, and the option to **download the RunInfoTable or Accession List**, in text format. The **RunInfoTable** is a very useful text summary of **all metadata for all runs** in the study, and the **Accession List** is a **list of all the SRR accession numbers for the study**.

Also on this page is a listing of each run and the corresponding sample it came from, as well as its associated metadata. This table is useful in that each row is "clickable", which allows you to select a subset of runs that you may be interested in. You'll notice that clicking a subset of runs spawns a new download option - a RunInfoTable & Accession List that is only relevant to your chosen subset.

<img src="../img/accession_list_red.png" width="600">

**Download the Accession list** for the data you are interested in to your desktop (everything included by default). 
**Copy the contents of this downloaded file to a new file on the cluster** using the following commands:

```bash
$ cd ~/Day5/SRA
```

During download, in addition to writing the fastq files, SRA-toolkit writes additional cache files, which are automatically directed to your home directory by default, even if you are working elsewhere. Because of this, we need to write a **short configuration file** to tell SRA-toolkit to **write its cache files to the scratch space**, instead of our home, to avoid running out of storage.

Given one single SRR, it is possible to convert that directly to a fastq file on the server, using [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/HowTo:-Access-SRA-Data) which is a toolkit created by NCBI. This should be already downloaded and installed on most clusters used for biomedical purposes.

```bash
$ fastq-dump SRR1013512
```

> **NOTE: Downloading Paired End Data:**
> Unlike the standard format for paired end data, where we normally find two fastq files labeled as `sample1_001.fastq` and `sample1_002.fastq`, **SRR files can be very misleading in that even paired end reads are found in one single file**, with sequence pairs concatenated alongside each other. Because of this format, **paired files need to be split down the middle** at the download step. 
>
> SRA toolkit has an option for this called `--split-files`. By using this, one single SRR file will download as `SRRxxx_1.fastq` and `SRRxxx_2.fastq`.
>
> Furthermore, there is a **helpful improvement** for this option called `--split-3`, which splits your SRR into 3 files: one for read 1, one for read 2, and one for any orphan reads (ie: reads that aren't present in both files). This is important for downstream analysis, as some aligners require your paired reads to be in sync (ie: present in each file at the same line number) and orphan reads can throw this order off.

The **second script loops through our list of SRRs**, and calls the first script from within the loop, passing it the next SRR in the list.

> **NOTE: SRRs from Multiple Studies:**
> Sometimes, in a publication, the relevant samples under study are given as sample numbers (GSM numbers), not SRRs, and sometimes belong to different GEO datasets (eg: different parts of a series, or separate studies for case and control experiments/data). If this is the case, download the RunInfoTables for each of the relevant studies as shown, selecting only the relevant GSMs/SRRs in the table before download, and copy them into one file. The starting point for the parallel fastq dump is a list of SRRs - so it does not matter if they came from different studies.

---
*This lesson has been modified by Shahin Shahsavari using materials from members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

