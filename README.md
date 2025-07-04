---
author: Bethany J. Allen,Joëlle Barido-Sottani
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

### TreeAnnotator

TreeAnnotator is used to summarise the posterior sample of trees to produce a maximum clade credibility tree and summarize the posterior estimates of other parameters that can be easily visualized on the tree (e.g. node height). This program is also useful for comparing a specific tree topology and branching times to the set of trees sampled in the MCMC analysis.

### R / RStudio

We will be using R to analyze the outputs of our analyses. RStudio provides a user-friendly graphical user interface to R that makes it easier to edit and run scripts. (It is not necessary to use RStudio for this tutorial.)

----

# Practical: Total Evidence Tutorial

In this tutorial we will infer a phylogeny for the Osmundaceae, or "Royal Fern", family of plants, using a combination of both morphological and genetic data. We will discuss the model elements and choices that are relevant, but recommend {% cite Mulvey2025a --file Total-Evidence-Tutorial/master-refs.bib %} as a recent review which explains these choices in more detail, and provides practical advice on compiling the necessary data.

The aim of this tutorial is to:
- Learn how to read morphological data into BEAST2;
- Understand the differences between how BEAST2 treats morphological and genetic data;
- Practice visualising total-evidence phylogenies.

## The data

The dataset we are using was published by {% cite Grimm2015 --file Total-Evidence-Tutorial/master-refs.bib %}. They used the dataset to investigate the impact of clock model choices on phylogenetic inference, between analyses conducted via "node dating" in BEAST2, "total-evidence dating" in MrBayes and "fossilised birth-death" dating in FDPPDiv. The genetic sequences comprise a concatenation of seven plastid DNA sites, collated by {% cite Metzgar2008 --file Total-Evidence-Tutorial/master-refs.bib %}. The morphological data describe characteristics of the plants' rhizomes, or underground stems, and are taken from {% cite Bomfleur2015 --file Total-Evidence-Tutorial/master-refs.bib %}.

The dataset describes 33 "operational taxonomic units", or biological entities which lie at the tips of our phylogeny. These include four species of _Leptopteris_, and three species of _Todea_, with the remaining tips belonging to the genus _Osmunda_. Most are independent species, while four are individuals which all belong to _Osmunda cinnamomea_: fossils from the Cretaceous of Canada, and the Neogene of Japan and the USA, alongside a genetic sequence for living _Osmunda cinnamomea_, each of which are represented as separate tips.

We will now explore the nuances of our two data input files: these are _Osmundaceae dna.nex_ for our genetic sequences, and _Osmundaceae morph.nex_ for our morphological data.

### Genetic sequences

First we will open the `nexus` file containing the genetic sequences in our text editor, to take a look at what the file contains. 

> Open `Osmundaceae_dna.nex` in your preferred text editor.

<figure>
	<a id="fig:1"></a>
	<img style="width:75%;" src="figures/dna_nexus.png" alt="">
	<figcaption>Figure 1: The nexus file containing the genetic     sequences.</figcaption>
</figure>

On the first line, we see a tag denoting that our file is of the `nexus` type. We then begin a `block` which is called `data`. The data dimensions are described, stating that `ntax = 33` and `nchar = 8628`, denoting that we have sequences with a length of 8628 base pairs for each of 33 OTUs, or tips. We are then given some information about the formatting of the sequences. The type of data is stated to be "dna". "Gaps" in the data are denoted using `-`, while "missing" values are denoted using `?`.

Following this header, the data itself is presented, following the tag `matrix`. Each tip label is stated, giving the name of the tip and a number (more on this later), followed by the genetic sequence associated with the OTU. The first is `Leptopteris_fraseri_0`. For this tip, we can see values of G, C, A and T corresponding to our DNA base reads, interspersed with blocks of question marks corresponding to sections of the DNA strand which could not be sequenced.

Scrolling down to the bottom, we see that `Osmunda_plumites_160` has a sequence which only consists of question marks. This is because _Osmunda plumites_ is one of our fossil tips. While we do not have a genetic sequence for this species, it must still be included in our nexus file containing the genetic sequences, but in a way which informs our model that we do not know the DNA bases that a sample from this species would contain.

<figure>
	<a id="fig:2"></a>
	<img style="width:75%;" src="figures/Osmunda_plumites.png" alt="">
	<figcaption>Figure 2: Extinct OTUs with no genetic data are represented with a sequence of question marks.</figcaption>
</figure>

To complete the file, the phrase `end;` is used to denote the end of the `data` block. No other blocks are needed. This example file shows one of the simplest ways of formatting a `nexus` file for input into BEAST2; if issues are encountered when reading your own sequence data into BEAST2, it is recommended to use the simplest possible formatting which still encodes the data.

### Morphological data

We will now look at the `nexus` file containing the morphological data for comparison.

> Open `Osmundaceae_morph.nex` in your preferred text editor.

<figure>
	<a id="fig:3"></a>
	<img style="width:75%;" src="figures/morph_nexus.png" alt="">
	<figcaption>Figure 3: The nexus file containing the morphological data.</figcaption>
</figure>

The first thing to note is that the overall format of the `nexus` file is very similar to that for our DNA sequences. We see a header, followed by a matrix. However, it is clear that this file is much shorter. For each of our DNA sequences we had a total of 8628 characters, but here we can see that `nchar = 25`. This means that our morphological data encodes 25 different traits possessed by each of our tips.

This time, the data format is described as "standard" rather than "dna". When using a "dna" format, we automatically determine that our sequences will consist of the letters G, C, A and T. But for this "standard" format, there is no clear character convention, and so we must describe which characters are contained within the matrix. This is done here using `symbols = "012"`. This means that the digits 0, 1 and 2 are used to denote the different character states, with no other characters present in the matrix. We therefore know that in this dataset, any single character can only have a maximum of three possible states. As with the DNA sequences, gaps and missing data are denoted using `-` and `?` respectively.

Beneath this, again, we see the name of each OTU followed by its character states for each trait. This time, no tip has only question marks, so we have at least some morphological data for every tip. But overall, there are still a high number of question marks spread across the dataset, with unknown states making up a larger proportion of the data compared to the DNA. For example, we can see that for `Todea_papuana_0`, only five of the characters are coded, with all other states unknown.

Another thing to note is that not all of our morphological "sequences" are the same length. This is because some contain bracketed pairs of values. For example, we see that `Osmunda_cinnamomea_0` contains two, the first being `{1 2}`, and the second being `{0 1}`. Here, we are denoting *ambiguity*. For these particular characters, we are telling the model that the value could be either of these two states. It is important to note that the model will interpret this information to mean that the state could be *either one* of these values, *not both*. If the morphological character in question was not coded with this information in mind, you risk violating the model of morphological evolution.

Here we have highlighted a handful of ways in which morphological data must be carefully formatted to be compatible with BEAST2. Different phylogenetic inference software use different models for morphological character evolution, and it is important to bear in mind the relationship between the morphological dataset design and the assumptions made by the model in your software of choice.

## Creating the Analysis Files with BEAUti

We will now start creating our BEAST2 model in BEAUti.

### Importing the data

First, we will load in our genetic data. This should be familiar to you if you have completed any of the other tutorials.

> Navigate to **File > Import Alignment**. 
> Select `Osmundaceae_dna.nex` and select **Open**.

We can now see this dataset listed in our **Partitions** tab. We can see that the dataset contains 33 **Taxa** and 8628 **Sites**, and has the **Data Type** "nucleotide": this corresponds to the information we saw in the `nexus` header, so it seems the data have been read in correctly.

<figure>
	<a id="fig:4"></a>
	<img style="width:75%;" src="figures/gen_sum.png" alt="">
	<figcaption>Figure 4: The summary of the genetic data displayed in BEAUti.</figcaption>
</figure>

Next, we will read in our morphological data.

> Navigate to **File > Add Morphological Data**. 
> Select `Osmundaceae_morph.nex` and select **Open**.

We now see a window asking whether we would like to "condition on recording variable characters only (Mkv)". This choice determines how the model treats the data with respect to **invariant sites**.

<figure>
	<a id="fig:5"></a>
	<img style="width:50%;" src="figures/variable_chars.png" alt="">
	<figcaption>Figure 5: The window asking whether to condition on recording only variable characters.</figcaption>
</figure>

What does this mean? We will start by considering our DNA sequences. Here, our data represent a complete sequence for each OTU, telling us which DNA base has been read at each position in the sequence (except for where we have missing values, but these are also included in the matrix, with a "?"). Imagine that for the first base in our sequences, all of the OTUs have the value "A". To use the technical term, this site is **invariant**. If all of the OTUs have the same value, we cannot use this particular position in the sequences to place the tips into two or more subgroups, and therefore to put our tips into clusters to inform the phylogenetic topology: it is **phylogenetically uninformative**. However, invariant sites can be very important with regards to **time**. We know that invariant sites have probably never changed state across any of our branches throughout the course of the evolutionary process we are modelling. When we infer a per-site rate of substitution across our sequences, invariant sites make up an important part of this calculation.

Now let's translate this thinking to morphological data. What do these data represent? Each character describes a specific trait, with the values given to each tip denoting the state in which that character exists in that OTU. The meaning of our morphological characters and states are much less rigid than for our genetic data {% cite Goloboff2019 --file Total-Evidence-Tutorial/master-refs.bib %}. The most relevant manifestation of this point here is that it is much easier to describe and code characters which differ between our OTUs, in comparison with characteristics that are shared by all of them. This means that morphological datasets are typically biased towards morphological characters which vary, with many not containing **any** invariant characters. This is known as **ascertainment bias**. As we mentioned for the genetic data, invariant characters form an important part of the calculation of per-character rates of change, but these are not captured in our dataset. Because of this, phylogenetic models have been adapted to correct for the non-sampling of invariant sites in morphological data. This is the 'variable', or **v**, part of the **MKv** model of morphological character evolution {% cite Lewis2001 --file Total-Evidence-Tutorial/master-refs.bib %}; we will discuss what **Mk** means later.

The original Osmundaceae dataset used by {% cite Grimm2015 --file Total-Evidence-Tutorial/master-refs.bib %} actually contained 33 morphological characters, with eight of these being **invariant**. However, BEAST2 currently will not accept morphological data containing invariant sites unless character descriptions are provided; if you attempt to read a morphological dataset containing invariant sites into BEAUti, you will get an error stating that this is not possible. For the purposes of this tutorial, we removed these invariant characters from the dataset for you. In your own dataset, you will need to either remove the invariant characters or provide character descriptions (see section **Advanced topics** for more information on how to do that).

To put this knowledge into action, we now know that our dataset does not include invariant sites, and so we would like to use the **Mkv** model to correct for ascertainment bias.

> Select **Yes** to conditioning on variable characters.

We can see that we now have three data partitions in the **Partitions** tab. Beneath our genetic data, our morphological dataset has been split into two. `Osmundaceae_morph2` contains 21 sites, and `Osmundaceae_morph3` contains 9 sites. You will notice that the sum of these two values is 30, higher than the 25 characters we saw in the raw morphological data: this is because five additional characters have been added to compensate for the ascertainment bias. Each partition has its own site model, and we will leave this setting to the default. We will discuss what this means when we set up the **Site model** tab.

<figure>
	<a id="fig:6"></a>
	<img style="width:75%;" src="figures/all_dat_sum.png" alt="">
	<figcaption>Figure 6: The summary of the genetic and morphological data displayed in BEAUti.</figcaption>
</figure>

Next we will take a look at the clock models. For now, we can see that we have one clock for the DNA sequences, named `Osmundaceae_dna`, and a second which is used for both of the morphological partitions, named `Osmundaceae_morph`. This means that we consider all of our morphological data to have evolved under the same clock, with a unified rate. This seems a reasonable assumption, so we will leave this setting to the default.

The following column defines the tree settings. We want to infer a single tree which is informed by both the genetic and morphological data, so we will link all three trees.

> In the **Tree** column, click on the name **Osmundaceae_dna**. 
> Type in **tree** and hit enter. 
> Using the drop-down menus, change the tree to **tree** for both of the morphological partitions.

<figure>
	<a id="fig:7"></a>
	<img style="width:75%;" src="figures/linked_trees.png" alt="">
	<figcaption>Figure 7: The summary of the genetic and morphological data displayed in BEAUti, with trees now linked.</figcaption>
</figure>

We are now ready to move onto our tip date settings.

### Setting initial tip ages

We know that our dataset contains both living and extinct species, which means that we have tips which were sampled through time. As a result, we need to give the model this information.

> Select the **Tip Dates** tab, and check the box to **Use tip dates**.

We now see a table which lists the tip names, along with **Date** and **Age/Height** values. At present these latter two columns are currently all zero values. These values correspond to the **initial** ages these tips will have: for some they will be sampled and inferred over the course of the MCMC, but we will come back to this point later.

<figure>
	<a id="fig:8"></a>
	<img style="width:75%;" src="figures/default_tips.png" alt="">
	<figcaption>Figure 8: The default settings for the tip dates.</figcaption>
</figure>

While it is possible to enter these age values manually, it is easier to specify them within our `nexus` file, which is what we have done here. Each tip label ends with an underscore followed by a number, and this number is the age, in millions of years, which we would like for the initial age of that tip. We will read these into the table.

> Ensure that the switch at the top for **Dates specified** is set to **numerically as... year**. 
> Change the second drop-down menu from **Since some time in the past** to **Before the present**. 
> Then select **Auto-configure**.
>
> Ensure that the switch is set to **use everything**. 
> Change the drop-down to **after last**, and ensure that the text box contains a single underscore, `_`. 
> Select **OK**.

<figure>
	<a id="fig:9"></a>
	<img style="width:75%;" src="figures/configuring_auto_tips.png" alt="">
	<figcaption>Figure 9: The window for setting the tip dates to autoconfigure from the tip labels.</figcaption>
</figure>

The **Date** and **Age/Height** values should now be repopulated using the numbers in the tip labels.

<figure>
	<a id="fig:10"></a>
	<img style="width:75%;" src="figures/configured_tips.png" alt="">
	<figcaption>Figure 10: The correctly configured tip dates.</figcaption>
</figure>

### Site models

We will now set up the site models. These models describe the evolutionary processes which we think govern the changing of states in both our molecular and morphological datasets, so it is particularly important to understand the choices we are making in this tab when setting up a total-evidence analysis.

> Select the **Site Model** tab.

We can see that our three partitions are listed in the table on the left, with the genetic data listed first. Here we will keep this model simple, but we recommend dipping into the other tutorials for more advice on how to set up site models for genetic sequence data.

> Set the **Gamma Category Count** to **4**, leaving the settings for the shape of the gamma distribution to the defaults. 
> Change the substitution model to **HKY**, leaving the initial **Kappa** value at **2.0**.
> Set the **Kappa** and the **Frequencies** parameters to **estimated**.

<figure>
	<a id="fig:11"></a>
	<img style="width:75%;" src="figures/dna_site_model.png" alt="">
	<figcaption>Figure 11: The settings for the site model for the genetic sequence data.</figcaption>
</figure>

> Select **Osmundaceae_morph2** in the left-hand table.

We can now begin thinking about our model(s) of morphological evolution. We can see that the default **Subst Model** is the **Lewis Mk**. At present this is the only model of morphological evolution available in BEAST2. (If you check the drop-down for this box, there is also the option of a **Mutation Death Model**. This is a **Dollo model**, which handles exclusively binary character states. It is configured specifically for word matrices used in language evolution, so we will not discuss this model further here.)

<figure>
	<a id="fig:12"></a>
	<img style="width:75%;" src="figures/morph_site_model.png" alt="">
	<figcaption>Figure 12: The settings for the site model for the morphological data.</figcaption>
</figure>

In brief, what is the Lewis Mk model? It is the **Markov _k_-states model described by {% cite Lewis2001 --file Total-Evidence-Tutorial/master-refs.bib %}. For those who are familiar with genetic substitution models, it is the morphological equivalent of the Jukes-Cantor 69 model. The model assumes that all character state transitions are just as likely as each other. For genetic sequences, this means that an A base is equally likely to transition to a G, C, or T.

In morphological evolution, this means several important things. Our Lewis MK model assumes that the transition from 0 to 1 is just as likely as the transition from 1 to 0 {% cite Lewis2001 --file Total-Evidence-Tutorial/master-refs.bib %}. If a character represents whether a complex evolutionary feature is present or not, our transitions therefore represent acquisition and loss. We might think it is much harder to lose such a character than to gain it: that would be a violation of this model. Further, the model assumes that, for example, the transition from 0 to 1 is just as likely as the transition from 0 to 2. This means that we assume character states have been coded as independent entities, but if these states are actually ordered (0 must transition to 1, then 1 must transition to 2), this is also a violation of the model. Understanding how the morphological characters have been coded, and whether any violations exist between these thought processes and the way in which the site model operates, is therefore fundamental to determining how well the site model will perform, and potentially how accurate our inferred phylogeny is {% cite Goloboff2019 Mulvey2025b --file Total-Evidence-Tutorial/master-refs.bib %}.

Now we can return to our earlier question: why was the morphological data split into two partitions? Under the Lewis Mk model, we know that for each morphological character, our transition matrix is dependent on the number of states, _k_, which that character can take. As a result, `Osmundaceae_morph2` contains all characters with two possible states, and `Osmundaceae_morph3` contains all characters with three possible states. In order to determine these partitions, BEAST2 automatically uses the number of character states present in the dataset for each character. This means that if a character is only represented by 0 or 1 in the dataset, BEAST2 assumes that there is no "secret" state 2.

We now need to check our morphological site model parameters. When we read in the dataset, we already told the model that there are no invariant sites, so we can leave the **Proportion invariant** blank.

We also need to choose a value for the **Gamma Category Count**. This value determines how many transition rate categories we have, allowing for differences in transition rates between different morphological characters in the matrix. BEAST2 uses the recommended implementation for using models that condition on only sampling variable sites but also permit for among-character variation in evolutionary rates, so we can use gamma categories in combination with the Mk model {% cite Capobianco2025 --file Total-Evidence-Tutorial/master-refs.bib %}. We just set this value to 4 for our molecular data. But as we mentioned before, the nature of morphological data means that we have much more variation in what different morphological characters actually describe than between different sites in a genetic sequence; this might lead us to think that we should allow for more different rate categories for our morphological data. However, we also have a much smaller total number of characters, giving us much less statistical power to actually infer these rate values compared to genetic data. To keep our analysis simple, we will leave the Gamma Category Count values for the morphological data at 0. This means we will calculate a single transition rate for each of the morphological partitions.

### Clock models

The next tab describes our clock models.

> Select the **Clock model** tab.

In the left-hand table, we see our two clocks, one for the genetic data and one for the morphological data. Clicking between them, we can see that for both, the default settings are for a **Strict Clock** with a **Mean clock rate** starting at 1.0, but which is **estimated** during the MCMC. We are dealing with a small phylogeny and a small number of characters, particularly for our morphological data, which are unlikely to provide us with enough information to fit a more complex clock model. We will therefore retain these defaults and use strict clocks.

<figure>
	<a id="fig:13"></a>
	<img style="width:75%;" src="figures/clock.png" alt="">
	<figcaption>Figure 13: The default settings for the molecular clock model.</figcaption>
</figure>

### Priors

Now it's time to set up our priors.

> Select the **Priors** tab.

Our first task is to choose our tree prior, which determines the type of evolutionary model we will use, and therefore the type of tree we will infer. We know that our phylogeny contains both extant and extinct tips, and because our datasets are small, we will keep things simple and assume constant evolutionary rates through time. This means that we need to use the **Fossilized Birth-Death Model**.

> In the drop-down box next to **Tree**, select **Fossilized Birth-Death Model**.
>
> Click the triangle next to **Tree** to view the additional tree prior options.

<figure>
	<a id="fig:14"></a>
	<img style="width:75%;" src="figures/FBD_options.png" alt="">
	<figcaption>Figure 14: The tree prior options for the fossilized birth-death model.</figcaption>
</figure>

Here we see many different options for our model. We can see that the FBD model will be parameterised with a diversification rate, turnover rate, and fossil sampling proportion, each of which will be estimated. We also have a value for **rho**, the proportion of extant tips which are included within the dataset, which is currently fixed at 1.0 (complete extant sampling). At time of writing, [World Flora Online](https://wfoplantlist.org/taxon/wfo-7000000433-2024-12?page=1) considers there to be 28 accepted living species in Osmundaceae. Our dataset contains 15 species dated to the present day (check the number with an age of "0" in the **Tip Dates** tab): we will therefore set the rho value to 15/28 = 0.536.

> Set the **rho** value to **0.536**.

We will leave the remaining tree parameter priors to their defaults.

The next two priors describe the clock rate parameters for the molecular and morphological data. We expect these values to be fairly low, so we will change both to exponential priors with a mean of 1.0.

> In the drop-down box next to **clockRate.c:Osmundaceae_dna**, select an **Exponential** distribution. By default, this should have a mean of 1.0; you can check this using the triangle on the left to see the extended options. 
> Repeat this for **clockRate.c:Osmundaceae_morph**.

<figure>
	<a id="fig:15"></a>
	<img style="width:75%;" src="figures/exp_clock_rate.png" alt="">
	<figcaption>Figure 15: Setting the molecular clock rate to have an exponential prior distribution.</figcaption>
</figure>

We also have three parameters for our molecular site model: **freqParameter**, **gammaShape** and **kappa**. We will leave these at their default values.

We can now move on to our FBD parameters. First, we will consider the origin time, which describes how old the root of our phylogeny should be.

> Click the triangle to the left of **OriginFBD** to view the default prior.

Here we can see that the default prior is a **uniform distribution** between 0 and infinity, with a starting value of 100. Practically, this means that the only constraint being placed on our phylogeny is that the root has to be older than the age of our oldest fossil tip - our **Tip Dates** tab tells us that the oldest starting tip age is 184 million years ago.

It is important for us to choose a more meaningful prior here, to ensure that the timescale of our phylogeny fits our prior knowledge. For example, we can check the [Paleobiology Database](https://paleobiodb.org/#/), a large open-access database of fossil occurrences. A search for the family **Osmundaceae** tells us that the oldest fossils in the database are from the [Carboniferous](https://paleobiodb.org/classic/basicTaxonInfo?taxon_no=54780), which corresponds to an age of around 360 to 300 million years ago. This is substantially older than the oldest fossil in our phylogeny. It is also older than the root age inferred in the original paper by {% cite Grimm2015 --file Total-Evidence-Tutorial/master-refs.bib %}, which corresponds to the Permian-Triassic boundary, 250 million years ago. To permit all of these options but exclude more extreme values, we will change the limits on our **uniform distribution** to correspond to the start of the Carboniferous, and the middle of the Triassic.

> Change the **Upper** limit of the uniform distribution to **360**, and the **Lower** limit to **230**. 
> Click the **initial =** button to view the starting values, and change the **Value** to **250**, then click **OK**.

<figure>
	<a id="fig:16"></a>
	<img style="width:75%;" src="figures/origin_prior.png" alt="">
	<figcaption>Figure 16: The settings for the prior distribution on the origin.</figcaption>
</figure>

Finally, we also have the three priors on our diversification and turnover rates, and fossil sampling proportion. We expect our diversification rate to be small, and our sampling proportion to be even smaller, so we will change these to exponential priors, while leaving the turnover parameter with the default uniform prior.

> In the drop-down box next to **diversificationRateFBD**, select an **Exponential** distribution. 
> Check that this distribution has a **mean** of **1.0**.

<figure>
	<a id="fig:17"></a>
	<img style="width:75%;" src="figures/diversification_prior.png" alt="">
	<figcaption>Figure 17: The settings for the prior distribution on the diversification rate.</figcaption>
</figure>

>In the drop-down box next to **samplingProportionFBD**, select an **Exponential** distribution. 
> Use the triangle on the left to open the extended options. 
> Give the distribution a **mean** of **0.1**. 
> Click the **initial =** button to view the starting values, and change the **Value** to **0.1**, then click **Ok**.

<figure>
	<a id="fig:18"></a>
	<img style="width:75%;" src="figures/sampling_prior.png" alt="">
	<figcaption>Figure 18: The settings for the prior distribution on the sampling proportion.</figcaption>
</figure>

### Adding tip age priors

As well as the priors which are already described in the **Priors** tab, we also need to set priors on the tip ages of our extinct OTUs. This will allow the ages of these tips to be sampled over the course of the MCMC. To do this, we need to give the model our uncertainty range for the age of each fossil tip. If you have previously completed the **Divergence Time Estimation** tutorial, this will be familiar to you, but we will walk this through with _Todea tidwellii_ as an example, which we consider to be between 140 and 129 million years old.

> Click the **+ Add Prior** button at the bottom of the **Priors** tab page. 
> When asked about the type of prior, select the **Sampled Ancestors MRCA Prior**.

<figure>
	<a id="fig:19"></a>
	<img style="width:50%;" src="figures/sa_mrca_prior.png" alt="">
	<figcaption>Figure 19: The window giving the option of adding a Sampled Ancestors MRCA Prior.</figcaption>
</figure>

> In the **Taxon set label** box, enter **tidwellii**. 
> In the left-hand box, find **Todea_tidwellii_135**, select it, and click the **>>** button to move it into the right-hand table. 
> Click **Ok**.

<figure>
	<a id="fig:20"></a>
	<img style="width:75%;" src="figures/add_tidwellii.png" alt="">
	<figcaption>Figure 20: The taxon set editor window for adding a tip prior for _Todea tidwellii_.</figcaption>
</figure>

>Using the drop-down box next to the newly-created prior, select a **Uniform** distribution. 
> Click the left-hand triangle to view the additional options. 
> Set the **Lower** value to **129** and the **Upper** value to **140**. 
> Check the box at the bottom for **Tipsonly**.

<figure>
	<a id="fig:21"></a>
	<img style="width:75%;" src="figures/modify_tidwellii.png" alt="">
	<figcaption>Figure 21: The settings for the tip prior on _Todea tidwellii_.</figcaption>
</figure>

Using this process, you can now enter the tip priors for each of the extinct tips, using the age uncertainties given in the box below. We have included the tips which have a fixed age for the sake of completeness, but it is not necessary to add a prior for these tips, as their initial ages are set to this date and if no prior is added, they will remain fixed at this age throughout the MCMC.

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
|---|---|

---

> Create a tip age prior for each extinct tip.

<figure>
	<a id="fig:22"></a>
	<img style="width:75%;" src="figures/tip_priors.png" alt="">
	<figcaption>Figure 22: The tip priors.</figcaption>
</figure>

Note that this is a time-consuming process, especially for large phylogenies. In order to automate this process, we have also provided the R script `scripts/add_age_uncertainty.R`, which can automatically add age uncertainty to a BEAST2 XML file. To use this script, first set up an XML file with your complete analysis setup, containing all the elements except for the tip age priors. Then, create an R table similar to the one shown above, containing a column `taxa` with the name of each fossil and columns `min_age` and `max_age` for the lower and upper bounds of each age interval. You can then run the script to produce an XML file containing the age uncertainty for each tip.

### Setting up MCMC

The last step is to set up our MCMC options.

>Select the **MCMC** tab.

The only setting we will change here is the chain length.

>Increase the chain length to **30000000** iterations.
>
>Leave all other settings at their default values and save the file as `Osmundaceae.xml`.

<figure>
	<a id="fig:23"></a>
	<img style="width:75%;" src="figures/mcmc_options.png" alt="">
	<figcaption>Figure 23: The MCMC settings.</figcaption>
</figure>

In one final step before we close BEAUti, we will also run an analysis which **samples from the prior**, to allow easier comparison between the shape of our priors and posteriors.

>Check the box to **Sample from prior**.
>
>Leave all other settings at their default values and save the file as `Osmundaceae_sfp.xml`.

We are now ready to run our analyses.

>Open BEAST2 and select `Osmundaceae.xml`. Hit **Run** to start the analysis.
>
>Once this is complete, reopen BEAST2 and select `Osmundaceae_sfp.xml`. Hit **Run** to start the analysis.

## Checking the logs
		
Once the BEAST2 analyses have finished running, we can take a cursory look at them using **Tracer**.
		
>Open **Tracer** and load both `Osmundaceae.log` and `Osmundaceae_sfp.log`.
		
First, we will check whether our empirical chain has converged. Once you have selected `Osmundaceae.log` in the top-left table, you will see alongside each of the parameter names are their `mean` and `ESS` values. `ESS` stands for **effective sample size**, and is a metric commonly used to determine whether a Bayesian analysis has converged. Values over 200 are typically taken to denote convergence; if any of the parameter values are below this, then the chain should be run for more iterations prior to analysis of the results.

It is likely that given our chain length, some but not all of the parameters will have converged. You can (and should) also confirm this visually: select any parameter which has an ESS over 200, then click the `Trace` button at the top of the window, and you should see the characteristic "caterpillar" of a well-mixed chain, but select any parameter with an ESS below this value, and the trace will appear more undulating. For scientific study we should ensure complete convergence of all parameters, but for our purposes here we will examine the results as if they were converged.

<figure>
	<a id="fig:24"></a>
	<img style="width:75%;" src="figures/tracer_chain.png" alt="">
	<figcaption>Figure 24: A well-mixed chain trace, here for the parameter **treeLength**.</figcaption>
</figure>

To investigate the differences between our prior and posterior distributions, we will compare the empirical chain to the one which we sampled from the prior.

>Select both `Osmundaceae.log` and `Osmundaceae_sfp.log` in the top-left table
> Switch to the **Marginal Density** tab on the right to view the distributions from both chains on the same plot. 
> The **Legend** drop-down at the bottom can be used to add a legend, indicating which distribution is from which chain.

Exploring the different parameters, we can see that for some, the prior and posterior are highly similar. For example, this is the case for most of the fossil tip heights given at the bottom of the table. This means that our prior is highly influential on the values we see in the posterior. However, for some other parameters, this is not the case, such as for parameters linked to the molecular data (e.g. `gammaShape`, `kappa`, `freqParameter`). This means that the shape of the posterior is being driven by the data rather than just the prior.

<figure>
	<a id="fig:25"></a>
	<img style="width:75%;" src="figures/gamma_shape.png" alt="">
	<figcaption>Figure 25: Comparison of the sampled-from-prior and posterior gamma shape estimates.</figcaption>
</figure>

What about for our morphological data? If we select `clockRate.morph`, we see a blank graph. At first this seems highly concerning, but actually we can work out why: if we select the two logs individually, we see that the prior distribution ranges between 0 and around 7, while the posterior distribution consists a highly concentrated spike at a value of around 0.001. Tracer is simply failing to plot these two very different distributions on the same set of axes - which might be unhelpful for now, but does tell us that our posterior distribution for this parameter is being heavily influenced by the data.

## Checking the trees

Next we will create a summary tree, to investigate the phylogenies we have inferred. We will do this using `TreeAnnotator`.

>Open `TreeAnnotator`, from within the BEAST2.7.x folder.
>
>Keep the default summary settings (10% burn in, posterior probability limit of 0, **Maximum clade credibility tree**) except for the node heights, for which we will use **Median heights**.
>
>For the input tree file, use **Choose file...** to navigate to **Osmundaceae-tree.trees**. Change the name of the output file to **Osmundaceae-MCC.tree**.
>
>Click **Run**.

This process should be rapid. You may see a warning which flags that the trees contain **sampled ancestors**. This is an integral part of the fossilised birth-death model, whereby some fossils are placed directly onto the branches of later descendants. We will visualise this shortly.

>Using a web browser, open [Icytree](https://icytree.org). Drag and drop **Osmundaceae-MCC.tree** into the browser screen.

You should now be able to see your maximum clade credibility tree visualised. We can change some of the settings to better see the features of trees containing fossils.

>Using the top menu, select **Style > Axis > Age** to add a time axis to the phylogeny.
>
>In the **Style** menu, uncheck the **Collapse zero-length edges** option. Now our sampled ancestors are visualised in the tree, on zero-length branches.

<figure>
	<a id="fig:26"></a>
	<img style="width:75%;" src="figures/MCC_tree.png" alt="">
	<figcaption>Figure 26: The visualised MCC tree.</figcaption>
</figure>

A quick count shows us that in our MCC tree, five of our fossils are sampled ancestors. We can consider this in more detail by looking at the _Osmunda cinnamomea_ clade. Remember that our dataset was set up to include four tips pertaining to this species, three of which were fossils and one the living (genetic) sample. We might expect that our three fossils should be inferred as sampled ancestors of our living branch. How did our analysis do? We can see that our _Osmunda cinnamomea_ fossil from the Cretaceous of Canada is indeed retrieved as a sampled ancestor of the living species. However, the other two fossils, from the Neogene of the USA and Japan, cluster together in a neighbouring cherry. _Osmunda precinnamomea_ branches off between our Cretaceous fossil and the rest of the species.

Our MCC tree represents the single tree which contains the clades which are most supported within our MCMC chain. But it's important to look at uncertainty around this topology, to get a sense of how certain it really is. First we will look at uncertainty in the timescale of our tree.

>Using the top menu, select **Style > Node height error bars > height_95%_HPD**.

We can now see error bars attached to each node and tip, showing us the 95% highest posterior density for the age of that node or tip. Mousing over any branch also brings up a table which shows us the numbers attached to that bar.

Another way to investigate uncertainty in the tree topology is to quantify the proportion of the posterior containing each subtree within the MCC tree.

>Using the top menu, select **Style > Node height error bars > None** to remove these bars.
>
>Now select **Style > Internal node text > posterior**.

<figure>
	<a id="fig:27"></a>
	<img style="width:75%;" src="figures/MCC_tree_uncertainty.png" alt="">
	<figcaption>Figure 27: The MCC tree annotated with the proportion of the posterior containing each subtree.</figcaption>
</figure>

We can see that these values vary drastically across the tree. For example, we can see that relationships within _Leptopteris_, and between _Leptopteris_ and the extant species of _Todea_, are very certain, with almost all being very close to 1. By contrast, those relationships we looked at earlier within _Osmunda cinnamomea_ are much less certain.

It is commonplace within the phylogenetics literature to only present the MCC topology for which node support is at or above 50%. Nodes with support beneath this value are deemed too poorly supported to be reasonably interpreted, and so instead are generally collapsed into polytomies.

# Advanced topic: specifying the number of states

As mentioned in the **Site model** section of the tutorial, BEAUti can automatically estimate the number of states from the alignment. Note that this is based on the number of different states present in the alignment and not on numerical values: for instance, if a character has states 0, 1 and 3, this character will be counted as having **3** states, not 4.

However, it is also possible to provide the number of possible states for each character in the alignment file by writing a **CHARSTATELABELS** block. This block needs to contain a description for each character (one character per line), and is of the form
```
X LABEL / STATES,
```
where X is the index of the character, LABEL is the character description, and STATES is a list of all possible states, separated by spaces. The number of states for each character is then calculated by BEAUti based on the provided STATES list.

As mentioned earlier, by default BEAUti does not allow invariant characters in morphological character matrices, because there is no way for the software to automatically calculate how many states an invariant character is supposed to have, and so it cannot specify a transition matrix for these characters. However, invariant characters can be included into the alignment if the number of possible states for these characters is specified. 

An example Nexus file with character descriptions (including some invariant characters) is provided in `Penguins_charstatelabels.nex`.

# Advanced topic: ordered characters

As mentioned in the **Site model** section of the tutorial, the Mk Lewis model assumes that transitions between all character states are equally likely for a given character. However, some characters are **ordered**, meaning that a lineage cannot directly transition from state **0** to state **2**, but needs to go through state **1** first. This can happen when a character is discretized from a continuous trait, for instance body size: in this situation it is very unlikely that a lineage can go from the "small" category to the "large" category without going first through the "medium" category.

So can we account for ordered characters in our analysis? This cannot be done through BEAUti, but will require manually editing the XML file. Note that manually edited XML files can generally not be loaded again into BEAUti. When specifying such an analysis, it is thus a good idea to specify as much of the analysis as possible in BEAUti, and leave the manual changes for last.

The first step is to specify which characters should be considered as ordered. If all characters with the same number of states (for instance all ternary characters) are ordered, then we can use the default partition created by BEAUti, otherwise we will need to specify our own partition which contains only the ordered characters.

The second step is to define a rate matrix for our ordered characters, which will have a value of **1.0** for **possible** transitions, and **0.0** for **impossible** transitions. Note that this matrix will also depend on the number of states for the character, but we show here an example for a character with 3 states.

{% eq \left( \begin{array}{ccc}
      -1.0 & 1.0 & 0.0 \\
      1.0 & -2.0 & 1.0 \\
       0.0 & 1.0 & -1.0
    \end{array} \right) %}

We then need to convert this matrix into a BEAST2 object, by simply listing all elements in order (by row, left to right). Note that the diagonal elements are fully defined by the rest of the matrix (since all rows must equal to 0) and so do not need to be written. This results in the following XML code, which can be placed after the alignment (between the `</data>` and `<run>` elements):
```xml
<parameter id="3ordered_ratematrix" dimension="6" name="rates">1.0 0.0 1.0 1.0 0.0 1.0 </parameter>
```

The last step is to specify which likelihood calculation should use our ordered transition matrix. For this, we should first find the likelihood associated with our ordered partition. If we are using the default partition, it should look like this (where "penguins" is the name of the alignment and "3" the number of states for ordered characters):
```xml
<distribution branchRateModel="@RelaxedClockModel.c:penguins" id="morphTreeLikelihood.penguins4" spec="TreeLikelihood" useAmbiguities="true" tree="@tree">
    <data id="penguins3" spec="FilteredAlignment" ascertained="true" data="@penguins" excludefrom="50" excludeto="53"
        filter="5,8,12,15,31-32,34,44,56,58,71-74,76,83,88,96-97,104,109,117-118,120,122,126-127,137,140,143,151,156,158,160,167,169,171,186,194,199-200,206,209-210,214,218,222,227,235,238,246-248">
        <userDataType id="morphDataType.penguins3" spec="beast.base.evolution.datatype.StandardData" ambiguities="01 12 45" nrOfStates="3"/>
    </data>
    <siteModel gammaCategoryCount="4" id="morphSiteModel.s:penguins3" shape="@gammaShape.s:penguins" spec="SiteModel">
        <parameter estimate="false" id="mutationRate.s:penguins3" name="mutationRate">1.0</parameter>
        <parameter estimate="false" id="proportionInvariant.s:penguins3" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>
        <substModel datatype="@morphDataType.penguins3" id="LewisMK.s:penguins3" spec="morphmodels.evolution.substitutionmodel.LewisMK"/>
    </siteModel>
</distribution>
```

We then need to replace the default Mk Lewis substitution model by one using the transition matrix we just defined, i.e. replace this element:
```xml
        <substModel datatype="@morphDataType.penguins3" id="LewisMK.s:penguins3" spec="morphmodels.evolution.substitutionmodel.LewisMK"/>
```
with this one:
```xml
        <substModel id="morph3.substmodel" spec="GeneralSubstitutionModel" rates="@3ordered_ratematrix">
            <frequencies id="morph3_freqs" frequencies="0.333 0.333 0.334" spec="Frequencies"/>
        </substModel>
```
Note that we need to specify the frequencies for all states as well. Here we have chosen to leave them all equal (just as in the Mk Lewis model).

An example XML file with ordered characters is shown in `Penguins_ordered.xml`.

----

# Useful Links

- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file Total-Evidence-Tutorial/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- BEAST 1 website and documentation: [http://beast.bio.ed.ac.uk](http://beast.bio.ed.ac.uk)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users) 

----

# Relevant References

{% bibliography --cited --file Total-Evidence-Tutorial/master-refs.bib %}

