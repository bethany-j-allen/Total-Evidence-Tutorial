---
author: Bethany J. Allen
level: Intermediate
title: Total Evidence Tutorial
subtitle: Inferring phylogenies using genetic sequences and morphological characters
beastversion: 2.7.7
tracerversion: 1.7.x
---

# Background

In macroevolution, phylogenies provide a hypothesised roadmap describing the relationships between species (or higher taxa) and the timing of their divergences. Typically we aim to infer the most complete and accurate phylogeny possible, to obtain a reliable picture of the tempo and mode of evolution in that clade. Genetic sequences have become relatively easy to obtain for living species, particularly through openly accessible databases such as [NCBI Genbank](https://www.ncbi.nlm.nih.gov/genbank/). While the fossil record can provide us with valuable insights into species in deeper time, we cannot obtain genetic sequences for these species, and to integrate them into phylogenies we must instead rely on their morphological characteristics. Phylogenetic approaches which combine both genetic and morphological data (and typically both extant and extinct species) are known as **total-evidence** methods {% cite Ronquist2012 Zhang2016 --file Total-Evidence-Tutorial/master-refs.bib %}. Here we will demonstrate how to infer a total-evidence phylogeny using the fossilized birth death (FBD) model in BEAST2.

----

# Programs used in this Exercise 

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

[BEAST2](http://www.beast2.org) is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees {% cite Bouckaert2014 --file Total-Evidence-Tutorial/master-refs.bib %}. This tutorial uses the BEAST v{{ page.beastversion }}.

### BEAUti - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs, and the interface will be the same, on all computing platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### Any programmer-friendly text editor

We will need to edit the XML files produced by BEAUti, for which we'll need a text editor. It's best to use one designed for programmers as these include nice features such as syntax highlighting, which makes the code more reader-friendly. [Sublime Text](https://www.sublimetext.com) is a good option which is available for MacOS, Windows and Linux.

### Tracer

[Tracer](http://beast.community/tracer) is used to summarise the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v{{ page.tracerversion }}.

### R / RStudio

We will be using R to analyze the outputs of our analyses. RStudio provides a user-friendly graphical user interface to R that makes it easier to edit and run scripts. (It is not necessary to use RStudio for this tutorial.)

----

# Practical: Total Evidence Tutorial

In this tutorial we will infer a phylogeny for the Osmundaceae, or "Royal Fern", family of plants, using a combination of both morphological and genetic data. 

The aim of this tutorial is to:
- Learn how to read morphological data into BEAST2;
- Understand the differences between how BEAST2 treats morphological and genetic data;
- Practise visualising total-evidence phylogenies.

## The data

The dataset we are using was published by {% cite Grimm2015 --file Total-Evidence-Tutorial/master-refs.bib %}. They used the dataset to investigate the impact of clock model choices on phylogenetic inference, between analyses conducted via "node dating" in BEAST2, "total-evidence dating" in MrBayes and "fossilised birth-death" dating in FDPPDiv.

The dataset describes 33 "operational taxonomic units", or biological entities which lie at the tips of our phylogeny. These include four species of _Leptopteris_, and three species of _Todea_, with the remaining tips belonging to the genus _Osmunda_. Most are independent species, while four are individuals which all belong to _Osmunda cinnamomea_: fossils from the Cretaceous of Canada, and the Neogene of Japan and the USA, alongside a genetic sequence for living _Osmunda cinnamomea_, each of which are represented as separate tips.

We will now explore the nuances of our two data input files: these are _Osmundaceae dna.nex_ for our genetic sequences, and _Osmundaceae morph.nex_ for our morphological data.

### Genetic sequences

First we will open the nexus file containing the genetic sequences in our text editor, to take a look at what the file contains. 

>Open _Osmundaceae dna.nex_ in **Sublime Text** (or your preferred text editor).

On the first line, we see a tag denoting that our file is of the `nexus` type. We then begin a `block` which is called `data`. The data dimensions are described, stating that `ntax=33` and `nchar=8628`, denoting that we have sequences with a length of 8628 base pairs for each of 33 OTUs, or tips. We are then given some information about the formatting of the sequences.

### Morphological data

## Creating the Analysis Files with BEAUti

### Install BEAST2 packages

<figure>
	<a id="fig:1"></a>
	<img style="width:75%;" src="figures/BEAUti packages.png" alt="">
	<figcaption>Figure 1: The package manager window in BEAUti.</figcaption>
</figure>

>Open the **BEAST2 Package Manager** by navigating to **File > Manage Packages**.
>
>Install the **BDSKY** and **feast** packages by selecting them and clicking the **Install/Upgrade** button one at a time.

### Importing the data

### Setting initial tip ages

### Priors

### Adding tip age priors


| Fossil tip | Age uncertainty |
| --- | --- |
| _Todea tidwelli_ | 129--140 |
| _Osmunda cinnamomea_ Neogene USA | 16 |
| _Osmunda cinnamomea_ Neogene Japan | 12--14 |
| _Osmunda cinnamomea_ Cretaceous Canada | 72--83 |
| _Osmunda precinnamomea_ | 56--66 |
| _Osmunda chengii_ | 145--160 |
| _Osmunda bromeliaefolia_ | 11--14 |
| _Osmunda arnoldii_ | 56--66 |
| _Osmunda dowkeri_ | 56--59 |
| _Osmunda ilianensis_ | 4--6 |
| _Osmunda oregonensis_ | 34--40 |
| _Osmunda pluma_ | 56--66 |
| _Osmunda shimokawaensis_ | 12--14 |
| _Osmunda wehrii_ | 16 |
| _Osmunda pulchella_ | 183--185 |
| _Osmunda liaoningensis_ | 153--165 |
| _Osmunda wangii_ | 153--165 |
| _Osmunda plumites_ | 153--165 |
---


### Setting up MCMC

## Examining the results

We can now examine the results of our **fossilised-birth-death** model. Again, we need to read in the relevant log file and trim off the first 10% as burn-in.

```R
# Navigate to Session > Set Working Directory > Choose Directory (on RStudio)
# or change file name to the full path to the log file
#(Use "dinosaur_BDSKY_final.log" if you used our pre-cooked XML)
fbd_file <- "dinosaur_BDSKY.log"

#Read in coalescent log and trim 10% burn-in
fbd <- read.table(fbd_file, header = T) %>% slice_tail(prop = 0.9)
```

----

# Useful Links

- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file Total-Evidence-Tutorial/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- BEAST 1 website and documentation: [http://beast.bio.ed.ac.uk](http://beast.bio.ed.ac.uk)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users) 

----

# Relevant References

{% bibliography --cited --file Total-Evidence-Tutorial/master-refs.bib %}

