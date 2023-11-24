---
layout: default
title: FAQs
nav_order: 7
description: "Frequently Asked Questions on how to use and interprete microbetag"
---

## What is sensitive and heterogeneous in FlashWeave? 



![four_cases](../assets/images/flashweave_cases.png)



heterogeneity of these cross-study data sets, such as variation in habitats, measurement conditions, and sequencing technology, can lead to confounding associations, typically not addressed by current methods

Sensitive modes (-S) of FlashWeave use full abundance information (”continuous”), while fast modes (-F) work on discretized abundances. In contrast to FlashWeave, FlashWeaveHE excludes samples in which one partner is absent (colored grey). However, it still includes absences of OTUs within the conditioning sets.





Meta variables (MVs) are by default not normalized for FlashWeave-S and FlashWeaveHE-S and should thus, if necessary, be provided in a sensible pre-normalized format by the user. For FlashWeave-F and FlashWeaveHE-F, continuous MVs are by default discretized into two bins separated by their median.


## What about computing resources ? 


| taxonomy | \# of seqIds | time (sec) | way microbetag was ran       | 
|:--------:|:------------:|:----------:|:----------------------------:|
| GTDB     | 480          |  300       | manta included; on the fly 


