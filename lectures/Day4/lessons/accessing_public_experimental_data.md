Approximate time: 45 minutes

## Learning Objectives:

* Understand the type of data that is accessible from Gene Expression Omnibus (GEO)
* Demonstrate how to navigate the GEO website and FTP server
* Use the command-line interface to copy over data from GEO

# Accessing public NGS sequencing data

All types of next-generation sequencing (NGS) analyses require access to public data, regardless of whether we are analyzing our own data or the data output from someone else's experiment. Reference data is available online as well as the experimental data from many published (and unpublished) studies. To access this data generally requires the basic knowledge of the **command line** and an understanding about the **associated tools and databases**.

To find and download NGS experimental data and associated reference data, we will explore several different repositories. For **accessing experimental data**, we will explore the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) and the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) repositories. For **finding reference data**, we will navigate the [Ensembl database](http://useast.ensembl.org/index.html), [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) and [Wormbase](https://wormbase.org). We will focus on these repositories; however, other databases are also useful for exploring and downloading experimental and reference data, including [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) and the [NCBI Genome](https://www.ncbi.nlm.nih.gov/genome/), and we encourage you to explore these more on your own.

# Gene Expression Omnibus

To find public experimental sequencing data, the NCBI's Gene Expression Omnibus (GEO) website is a useful place to search. The requirement for many grants is that experimental data be uploaded to GEO and the sequence read archive (SRA); therefore, there are quite a few datasets on GEO available to search. The interface for finding data on GEO is relatively user-friendly and easily searchable.

## Finding and accessing data on GEO

### Searching GEO

To search GEO for particular types of data is relatively straight forward. Once on the [GEO website](https://www.ncbi.nlm.nih.gov/geo/) there are multiple different options for searching datasets. 

<img src="../img/geo_web.png" width="700">

The most straight-forward method can be found by clicking on 'Datasets' under the 'Browse Content' column. 

<img src="../img/geo_dataset.png" width="400">

The 'Datasets' link will open the GEO Dataset Browser; click on 'Advanced Search'.

<img src="../img/geo_browser.png" width="600">

All results will appear in a new window with clickable filters on the left-hand side. You can choose the filters, such as 'Organism' (human, mouse), 'Study type' (Expression profiling by high throughput sequencing), 'Publication dates' (1 year), etc. to filter the data for the desired attributes.

<img src="../img/geo_filter.png" width="700">

### Finding GEO data for a particular publication

To find data from a published paper on GEO, the paper will often provide the GEO accession number. For example, let's find the data associated with the paper, "MOV10 and FRMP regulate AGO2 association with microRNA recognition elements". First, we can navigate to the [article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4268400/).

<img src="../img/mov10_paper.png" width="600">

Then, we can search for the term **"GEO"**; different papers have different requirements for where this information is located. In this article, it is available in a separate section entitled "Accession Numbers".

<img src="../img/mov10_accession.png" width="600">

By clicking on the GEO accession number for the experiment of interest, the GEO page for this experiment will open.

<img src="../img/mov10_geo.png" width="700">

The GEO page contains information about the experiment, including:
	
- an experimental summary
- literature citation
- contact information
- links to the sample GEO pages
- link to the SRA project containing the raw FASTQ files
- metadata (matrix file)
- supplementary data
	
If we were interested in **downloading the raw counts matrix (`GSE50499_GEO_Ceman_counts.txt.gz`)**, which gives the number of reads/sequences aligning to each gene, then we could scroll down to **supplementary data** at the bottom of the page. 

We could download this file by clicking on the `ftp` link. In addition to the counts matrix file, we would probably also want the metadata for the file to know which sample belongs to which conditions by clicking on the "Series Matrix File(s)" link. 

> **NOTE:** The "Series Matrix" metadata file is a bit congested, and it may be easier accessing the metadata from the SRA instead as discussed later.

<img src="../img/mov10_download.png" width="600">

Now that we have these files, if we wanted to perform differential expression analysis, we could bring them into R to perform some data wrangling and analysis.

## Using the command line to download from GEO

While navigating the GEO website is a perfectly fine way to download data from GEO, oftentimes you may want to download multiple files at the same time or download large files that you don't want to store on your personal computer. Generally, the best way to do this is by downloading the files on the cloud.

### Downloading on a remote server 

We will now use the secure shell (`ssh`) command using the 'Terminal' program for Macs or 'Power Shell' for Windows (Use PuTTY if you don't have the latest updates). In the console from these programs, type:

```bash
$ ssh username@3.237.4.148
```

You will be prompted to enter your password.

This is the common way to access most servers in a secure manner.

Now, we can download our data to an appropriate directory. Good data management practices will ensure we have an organized project directory for our analysis. We can create and change directories to the folder to which we plan to download the data.

```bash
$ mkdir -p ~/access_public_data/GEO/counts

$ cd ~/access_public_data/GEO/counts
```

To access the FTP site, return to the [GEO home page](https://www.ncbi.nlm.nih.gov/geo/) and under the "Tools" header, click on "FTP site".

<img src="../img/geo_ftp.png" width="350">

This will take you to the directory to access all GEO data.

<img src="../img/geo_dir.png" width="350">

To download the data associated with the paper, "MOV10 and FMRP Regulate AGO2 Association with MicroRNA Recognition Elements", use the GEO ID given in the paper, `GSE50499`.

1. Navigate the FTP site to the `series/` folder
2. Find the `GSE50nnn/` directory 
3. Enter the `GSE50499/` folder
4. Use the `wget` command followed by the link to the `suppl/` directory (right-clicking and choosing 'Copy Link Address'). 

	<img src="../img/geo_folder_cp.png" width="500">

	Using the `wget` command to copy this dataset.

	```bash
	$ wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50499/suppl/GSE50499_GEO_Ceman_counts.txt.gz 
	```

5. You could also load the entire directory using:

	```bash
	$ wget --recursive --no-parent -nd ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50499/suppl/
	```
---

**Exercises**

1. What command would we use to download the data from this [paper](https://immunology.sciencemag.org/content/5/49/eabd1554.long): 

	<img src="../img/study_exercise_head.png" width="500">
	
	>
	>_**NOTE:** If you wanted to unpack the `.tar` file and decompress the `.gz` files, you could use the following commands:_
	>
	>```bash
	>$ tar -xvf GSE111889_RAW.tar 
	>
	>$ for file in *.gz; do gunzip $all; done
	>```

---

### Downloading on a local computer

If we are downloading a small file(s) to use on our personal computer, then it makes sense to directly download to our computer. Unfortunately, `wget` is not automatically available using the Mac Terminal or using GitBash on Windows. However, the `curl` command is available to transfer data from (or to) a server. We can  **download individual files** using `curl` by connecting to the FTP:

```bash
curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50499/suppl/GSE50499_GEO_Ceman_counts.txt.gz
```

Also, it's worth noting that we don't need to navigate the FTP site to find individual files to download, since the link on the GEO site should list a link to the file. By right-clicking on the `ftp` link on GEO, you can copy the 'ftp address' to use with the `wget` or `curl` command.

<img src="../img/geo_ftp_cl.png" width="700">

---
*This lesson has been developed by Shahin Shahsavvari using resources from members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
