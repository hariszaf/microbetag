---
layout: default
title: Run from Python
nav_order: 4
---


# Run microbetag from Python 


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



Similarly, in Python:

```python
>>> import requests
>>> url = "https://msysbio.gbiomed.kuleuven.be/ncbiTaxId-to-genomeId/853"
>>> r = requests.get(url)
>>> r.status_code
200
>>> r.json()
{"537007": ["GCA_002222595.2"]}
```



Convert your abundance table in a JSON object 

Make a dictionary with your arguments. 

```python

import requests
import json

data = json.load(open(), "r")

```




