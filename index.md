---
layout: default
title: Home
nav_order: 1
description: " A place with the background and the *how to* of the *microbetag* tool"
permalink: /
---

# annotating microbial co-occurrence networks
{: .fs-9 }

A place with the background and the *how to* of the *microbetag* tool
{: .fs-6 .fw-300 }

Under development
{: .label .label-yellow }



![microbetag logo](/assets/images/microbetag_logo.png){: width=15% }


<!--The ': .btn' flag denotes the button 
# The ': .fs-5' flag denotes the font size
# The ': .mb-4' flag denotes the margin-bottom: https://pmarsceill.github.io/just-the-docs/docs/utilities/layout/#spacing
# the 'mb'is the margin-bottom as said,  while the 'md' stands for a [responsive modifier](https://pmarsceill.github.io/just-the-docs/docs/utilities/responsive-modifiers/#responsive-modifiers)-->

[View it on GitHub](https://github.com/hariszaf/microbetag){: .btn .btn-purple .fs-5 .mb-4 .mb-md-0 }
[Image on DockerHub](https://hub.docker.com/r/hariszaf/microbetag){: .btn .btn-green .fs-5 .mb-4 .mb-md-0 }




---

## About

Microbial interactions play a fundamental role in deciphering the underlying mechanisms that govern ecosystem functioning. 
Co-occurrence networks have been widely used for inferring microbial associations or/and interactions from metagenomic data. 
However, spurious associations and tool - dependence confine the network inference. 
The integration of previous evidence or/and knowledge can increase or decrease the confidence level of the retrieved associations. 
This way, associations can be further investigated and more reliable conclusions can be drawn.  


*microbetag* is a command line tool that implements data integration techniques to enhance microbial co-occurrence network analysis for amplicon data. 
*microbetag* retrieves the KEGG modules that have been assigned to each of the species found related. 
Based on the pathway complementarity concept, pathways found in both taxa of an association are further explored to check whether the processes of each of the two taxa are complementary denoting a  positive interaction. 
Likewise, if the same processes are found to occur in both taxa, a negative interaction will be derived.

On top of that, *microbetag* integrates phenotypic information thanks to resources such as [FAPROTAX](https://github.com/knights-lab/BugBase), [BugBase](https://github.com/knights-lab/BugBase) and [Picrust2](https://github.com/picrust/picrust2); 
a series of environmental variables (pH optima, oxygen tolerance etc.) are assembled in each node of the network.
Their comparison in each pair of correlated taxa evaluates their corresponding association further. 




## Dependencies
To run *microbetag* you need to have [Docker](https://www.docker.com/) on your computing environment. 
As described from IBM, Docker is an open source containerization platform. 
It enables developers to package applications into containers—standardized executable components combining application source code with the operating system libraries and dependencies required to run that code in any environment.

You can install Docker in Linux, MaxOS or Windows systems by following the instructions you will finde [here](https://docs.docker.com/get-docker/).


### Get

Once Docker is available, to get *microbetag* you need to *pull* it from DockerHub. 
To do this, you need to run: 

```bash=
docker push hariszaf/microbetag
```

This way, the latest version of *microbetag* will be pulled. 
You may specify which version of *microbetag* you wish to pull by running instead:

```bash=
docker push hariszaf/microbetag:tagname
```
where `tagname` is the name of the specific version. 





## Run 

To run *microbetag* you need an OTU-table which follows one of the following sets of features: 

- taxonomy assignment using Silva 132 
- 



## About the project



## License

- *microbetag* is under GNU General Public License v3.0. For third-party components separate licenses apply. 

- This web-page is based on the Just the Docs theme that is distributed by an [MIT license](https://github.com/pmarsceill/just-the-docs/tree/master/LICENSE.txt).

