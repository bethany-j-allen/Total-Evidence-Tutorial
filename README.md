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

First we will open the `nexus` file containing the genetic sequences in our text editor, to take a look at what the file contains. 

>Open `Osmundaceae_dna.nex` in your preferred text editor.

On the first line, we see a tag denoting that our file is of the `nexus` type. We then begin a `block` which is called `data`. The data dimensions are described, stating that `ntax = 33` and `nchar = 8628`, denoting that we have sequences with a length of 8628 base pairs for each of 33 OTUs, or tips. We are then given some information about the formatting of the sequences. The type of data is stated to be "dna". "Gaps" in the data are denoted using `-`, while "missing" values are denoted using `?`.

Following this header, the data itself is presented, following the tag `matrix`. Each tip label is stated, giving the name of the tip and a number (more on this later), followed by the genetic sequence associated with the OTU. The first is `Leptopteris_fraseri_0`. For this tip, we can see values of G, C, A and T corresponding to our DNA base reads, interspersed with blocks of question marks corresponding to sections of the DNA strand which could not be sequenced.

Scrolling down to the bottom, we see that `Osmunda_plumites_160` has a sequence which only consists of question marks. This is because _Osmunda plumites_ is one of our fossil tips. While we do not have a genetic sequence for this species, it must still be included in our nexus file containing the genetic sequences, but in a way which informs our model that we do not know the DNA bases that a sample from this species would contain.

To complete the file, the phrase `end;` is used to denote the end of the `data` block. No other blocks are needed. This example file shows one of the simplest ways of formatting a `nexus` file for input into BEAST2; if issues are encountered when reading your own sequence data into BEAST2, it is recommended to use the simplest possible formatting which still encodes the data.

### Morphological data

We will now look at the `nexus` file containing the morphological data for comparison.

>Open `Osmundaceae_morph.nex` in your preferred text editor.

The first thing to note is that the overall format of the `nexus` file is very similar to that for our DNA sequences. We see a header, followed by a matrix. However, it is clear that this file is much shorter. For each our DNA sequences we had a total of 8628 characters, but here we can see that `nchar = 25`. This means that our morphological data encodes 25 different traits possessed by each of our tips.

This time, the data format is described as "standard" rather than "dna". When using a "dna" format, we automatically determine that our sequences will consist of the letters G, C, A and T. But for this "standard" format, there is no clear character convention, and so we must describe which characters are contained within the matrix. This is done here using `symbols = "012"`. This means that the digits 0, 1 and 2 are used to denote the different character states, with no other characters present in the matrix. We therefore know that in this dataset, any single character can only have a maximum of three possible states. As with the DNA sequences, gaps and missing data are denoted using `-` and `?` respectively.

Beneath this, again, we see the name of each OTU followed by its character states for each trait. This time, no tip has only question marks, so we have at least some morphological data for every tip. But overall, there are still a high number of question marks spread across the dataset, with unknown states making up a larger proportion of the data compared to the DNA. For example, we can see that for `Todea_papuana_0`, only five of the characters are coded, with all other states unknown.

Another thing to note is that not all of our morphological "sequences" are the same length. This is because some contain bracketed pairs of values. For example, we see that `Osmunda_cinnamomea_0` contains two, the first being `{1 2}`, and the second being `{0 1}`. Here, we are denoting *ambiguity*. For these particular characters, we are telling the model that the value could be either of these two states. It is important to note that the model will interpret this information to mean that the state could be *either one* of these values, *not both*. If the morphological character in question was not coded with this information in mind, you risk violating the model of morphological evolution.

The original dataset used by {% cite Grimm2015 --file Total-Evidence-Tutorial/master-refs.bib %} contained 33 morphological characters. However, eight of these were *invariant*, meaning that all of the OTUs had the same character state. BEAST2 currently will not accept morphological data containing invariant sites; if you attempt to read a morphological dataset containing invariant sites into BEAUti, you will get an error stating that this is not possible. For the purposes of this tutorial, we removed these invariant characters from the dataset, but this is something you may need to conduct yourself if attempting to read a dataset into BEAUti which was not originally set up for BEAST2.

Here we have highlighted a handful of ways in which morphological data must be carefully formatted to be compatible with BEAST2. Different phylogenetic inference software use different models for morphological character evolution, and it is important to bear in mind the relationship between the morphological dataset design and the assumptions made by the model in your software of choice.

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

