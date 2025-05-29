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

The first thing to note is that the overall format of the `nexus` file is very similar to that for our DNA sequences. We see a header, followed by a matrix. However, it is clear that this file is much shorter. For each of our DNA sequences we had a total of 8628 characters, but here we can see that `nchar = 25`. This means that our morphological data encodes 25 different traits possessed by each of our tips.

This time, the data format is described as "standard" rather than "dna". When using a "dna" format, we automatically determine that our sequences will consist of the letters G, C, A and T. But for this "standard" format, there is no clear character convention, and so we must describe which characters are contained within the matrix. This is done here using `symbols = "012"`. This means that the digits 0, 1 and 2 are used to denote the different character states, with no other characters present in the matrix. We therefore know that in this dataset, any single character can only have a maximum of three possible states. As with the DNA sequences, gaps and missing data are denoted using `-` and `?` respectively.

Beneath this, again, we see the name of each OTU followed by its character states for each trait. This time, no tip has only question marks, so we have at least some morphological data for every tip. But overall, there are still a high number of question marks spread across the dataset, with unknown states making up a larger proportion of the data compared to the DNA. For example, we can see that for `Todea_papuana_0`, only five of the characters are coded, with all other states unknown.

Another thing to note is that not all of our morphological "sequences" are the same length. This is because some contain bracketed pairs of values. For example, we see that `Osmunda_cinnamomea_0` contains two, the first being `{1 2}`, and the second being `{0 1}`. Here, we are denoting *ambiguity*. For these particular characters, we are telling the model that the value could be either of these two states. It is important to note that the model will interpret this information to mean that the state could be *either one* of these values, *not both*. If the morphological character in question was not coded with this information in mind, you risk violating the model of morphological evolution.

Here we have highlighted a handful of ways in which morphological data must be carefully formatted to be compatible with BEAST2. Different phylogenetic inference software use different models for morphological character evolution, and it is important to bear in mind the relationship between the morphological dataset design and the assumptions made by the model in your software of choice.

## Creating the Analysis Files with BEAUti

We will now start creating our BEAST2 model in BEAUti.

### Install BEAST2 packages

The first step to creating our model is ensuring that we have all of the BEAST2 packages we need.

<figure>
	<a id="fig:1"></a>
	<img style="width:75%;" src="figures/BEAUti packages.png" alt="">
	<figcaption>Figure 1: The package manager window in BEAUti.</figcaption>
</figure>

>Open the **BEAST2 Package Manager** by navigating to **File > Manage Packages**.
>
>Install the **sa** packages by selecting it and clicking the **Install/Upgrade** button.
>
>Close and restart **BEAUti**.

### Importing the data

First, we will load in our genetic data. This should be familiar to you if you have completed any of the other tutorials.

>Navigate to **File > Import Alignment**. Select `Osmundaceae_dna.nex` and select **Open**.

We can now see this dataset listed in our **Partitions** tab. We can see that the dataset contains 33 **Taxa** and 8628 **Sites**, and has the **Data Type** "nucleotide": this corresponds to the information we saw in the `nexus` header, so it seems the data have been read in correctly.

Next, we will read in our morphological data.

>Navigate to **File > Add Morphological Data**. Select `Osmundaceae_morph.nex` and select **Open**.

We now see a window asking whether we would like to "condition on recording variable characters only (Mkv)". This choice determines how the model treats the data with respect to **invariant sites**.

What does this mean? We will start by considering our DNA sequences. Here, our data represent a complete sequence for each OTU, telling us which DNA base has been read at each position in the sequence (except for where we have missing values, but these are also included in the matrix, with a "?"). Imagine that for the first base in our sequences, all of the OTUs have the value "A". To use the technical term, this site is **invariant**. If all of the OTUs have the same value, we cannot use this particular position in the sequences to place the tips into two or more subgroups, and therefore to put our tips into clusters to inform the phylogenetic topology: it is **phylogenetically uninformative**. However, invariant sites can be very important with regards to **time**. We know that invariant sites have probably never changed state across any of our branches throughout the course of the evolutionary process we are modelling. When we infer a per-site rate of substitution across our sequences, invariant sites make up an important part of this calculation.

Now let's translate this thinking to morphological data. What do these data represent? Each character describes a specific trait, with the values given to each tip denoting the state in which that character exists in that OTU. The meaning of our morphological characters and states are much less rigid than for our genetic data. The most relevant manifestation of this point here is that it is much easier to describe and code characters which differ between our OTUs, in comparison with characteristics that are shared by all of them. This means that morphological datasets are typically biased towards morphological characters which vary, with many not containing **any** invariant characters. This is known as **ascertainment bias**. As we mentioned for the genetic data, invariant characters form an important part of the calculation of per-character rates of change, but these are not captured in our dataset. Because of this, phylogenetic models have been adapted to correct for the non-sampling of invariant sites in morphological data. This is termed the **MKv** model of morphological character evolution.

The original Osmundaceae dataset used by {% cite Grimm2015 --file Total-Evidence-Tutorial/master-refs.bib %} actually contained 33 morphological characters, with eight of these being **invariant**. However, because of the rarity of having this information in datasets, BEAST2 currently will not accept morphological data containing invariant sites; if you attempt to read a morphological dataset containing invariant sites into BEAUti, you will get an error stating that this is not possible. For the purposes of this tutorial, we removed these invariant characters from the dataset for you, but this is something you may need to conduct yourself if attempting to read a dataset into BEAUti which was not originally set up for BEAST2.

To put this knowledge into action, we now know that our dataset does not include invariant sites, and so we would like to use the **Mkv** model to correct for ascertainment bias.

>Select **Yes** to conditioning on variable characters.

We can see that we now have three data partitions in the **Partitions** tab. Beneath our genetic data, our morphological dataset has been split into two. `Osmundaceae_morph2` contains 21 sites, and `Osmundaceae_morph3` contains 9 sites. Each partition has its own site model, and we will leave this setting to the default. We will discuss what this means when we set up the **Site model** tab.

Next we will take a look at the clock models. For now, we can see that we have one clock for the DNA sequences, named `Osmundaceae_dna`, and a second which is used for both of the morphological partitions, named `Osmundaceae_morph`. This means that we consider all of our morphological data to have evolved under the same clock, with a unified rate. This seems a reasonable assumption, so we will leave this setting to the default.

The last thing we need to look at in this tab is the tree settings. We want to infer a single tree which is informed by both the genetic and morphological data, so we will link all three trees.

>In the **Tree** column, click on the name **Osmundaceae_dna**. Type in **tree** and hit enter. Using the drop-down menus, change the tree to **tree** for both of the morphological partitions.

We are now ready to move onto our tip date settings.

### Setting initial tip ages

We know that our dataset contains both living and extinct species, which means that we have tips which were sampled through time. As a result, we need to give the model this information.

>Select the **Tip Dates** tab, and check the box to **Use tip dates**.

We now see a table which lists the tip names, along with **Date** and **Age/Height** values. At present these latter two columns are currently all zero values. These values correspond to the **initial** ages these tips will have: for some they will be sampled and inferred over the course of the MCMC, but we will come back to this point later.

While it is possible to enter these age values manually, it is easier to specify them within our `nexus` file, which is what we have done here. Each tip label ends with an underscore followed by a number, and this number is the age, in millions of years, which we would like for the initial age of that tip. We will read these into the table.

>Ensure that the switch at the top for **Dates specified** is set to **numerically as... year**. Change the second drop-down menu from **Since some time in the past** to **Before the present**. Then select **Auto-configure**.
>
>Ensure that the switch is set to **use everything**. Change the drop-down to **after last**, and ensure that the text box contains a single underscore, `_`. Select **OK**.

The **Date** and **Age/Height** values should now be repopulated using the numbers in the tip labels.

### Site models

We will now set up the site models. These models describe the evolutionary processes which we think govern the changing of states in both our molecular and morphological datasets, so it is particularly important to understand the choices we are making in this tab when setting up a total-evidence analysis.

>Select the **Site Model** tab.

We can see that our three partitions are listed in the table on the left, with the genetic data listed first. Here we will keep this model simple, but we recommend dipping into the other tutorials for more advice on how to set up site models for genetic sequence data.

>Set the **Gamma Category Count** to **4**, leaving the settings for the shape of the gamma distribution to the defaults. Check the box to **estimate** the **Proportion invariant**. Change the substitution model to **HKY**, leaving the initial **Kappa** value at **2.0** but **estimated**, as well as the **Frequencies**.
>
>Select **Osmundaceae_morph2** in the left-hand table.

We can now begin thinking about our model(s) of morphological evolution. We can see that the default **Subst Model** is the **Lewis MK**. At present this is the only model of morphological evolution available in BEAST2. (If you check the drop-down for this box, there is also the option of a **Mutation Death Model**. This is a **Dollo model**, which handles exclusively binary character states. It is configured specifically for word matrices used in language evolution, so we will not discuss this model further here.)

In brief, what is the Lewis MK model? For those who are familiar with genetic substitution models, it is the morphological equivalent of the Jukes-Cantor 69 model. The model assumes that all character state transitions are just as likely as each other. For genetic sequences, this means that an A base is equally likely to transition to a G, C, or T.

In morphological evolution, this means several important things. Our Lewis MK model assumes that the transition from 0 to 1 is just as likely as the transition from 1 to 0. If a character represents whether a complex evolutionary feature is present or not, our transitions therefore represent acquisition and loss. We might think it is much harder to lose such a character than to gain it: that would be a violation of this model. Further, the model assumes that, for example, the transition from 0 to 1 is just as likely as the transition from 0 to 2. This means that we assume character states have been coded as independent entities, but if these states are actually ordered (0 must transition to 1, then 1 must transition to 2), this is also a violation of the model. Understanding how the morphological characters have been coded, and whether any violations exist between these thought processes and the way in which the site model operates, is therefore fundamental to determining how well the site model will perform, and potentially how accurate our inferred phylogeny is.

Now we can return to our earlier question: why was the morphological data split into two partitions? Under the Lewis MK model, we know that for each morphological character, our transition matrix is dependent on the number of states which that character can take. As a result, `Osmundaceae_morph2` contains all characters with two possible states, and `Osmundaceae_morph3` contains all characters with three possible states. In order to determine these partitions, BEAST2 automatically uses the number of character states present in the dataset for each character. This means that if a character is only represented by 0 or 1 in the dataset, BEAST2 assumes that there is no "secret" state 2.

We now need to set up our morphological site model parameters. When we read in the dataset, we already told the model that there are no invariant sites, so we can leave the **Proportion invariant** blank. Because this dataset is so small, we will set our gamma category count to 2.

>Set the **Gamma Category Count** to **2**, leaving the settings for the shape of the gamma distribution to the defaults.
>
>Select **Osmundaceae_morph3** in the left-hand table, again setting the **Gamma Category Count** to **2**.

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

