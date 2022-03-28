---
layout: default
title: Home
nav_order: 1
description: " A place with the background and the *how to* of the *microbetag* tool - currently in production"
permalink: /
---

# annotating microbial co-occurrence networks
{: .fs-9 }

A place with the background and the *how to* of the *microbetag* tool - currently in production
{: .fs-6 .fw-300 }

<!--The ': .btn' flag denotes the button 
# The ': .fs-5' flag denotes the font size
# The ': .mb-4' flag denotes the margin-bottom: https://pmarsceill.github.io/just-the-docs/docs/utilities/layout/#spacing
# the 'mb'is the margin-bottom as said,  while the 'md' stands for a [responsive modifier](https://pmarsceill.github.io/just-the-docs/docs/utilities/responsive-modifiers/#responsive-modifiers)-->

[View it on GitHub](https://github.com/hariszaf/microbetag){: .btn .btn-purple .fs-5 .mb-4 .mb-md-0 }

[Image on DockerHUb](https://hub.docker.com/r/hariszaf/microbetag){: .btn .btn-green-100 .fs-5 .mb-4 .mb-md-0 }

---

##  Getting started

### Dependencies
To run *microbetag* you need to have [Docker](https://www.docker.com/) on your computing environment. 
As described from IBM, Docker is an open source containerization platform. 
It enables developers to package applications into containers—standardized executable components combining application source code with the operating system libraries and dependencies required to run that code in any environment.

You can install Docker in Linux, MaxOS or Windows systems by following the instructions you will finde [here](https://docs.docker.com/get-docker/).


### Installation

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


## About the project



## License

- *microbetag* 

- This web-page is based on the Just the Docs theme that is distributed by an [MIT license](https://github.com/pmarsceill/just-the-docs/tree/master/LICENSE.txt).

