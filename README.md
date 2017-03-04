# imGLAD

imGLAD is a computational tool for detection of bacterial genomes in metagenomic datasets. For license information, see
[LICENSE](./LICENSE).

## Overview

The software consists of two parts the first part creates a series of metagenomic datasets, the datasets are created in such a way that the target organism is present in half of them and absent in the other half.

##System requirements

Python 2.7 or higher.

ART 2.5.8 or higher.

Either BLAST 2.2.28 (or higher) or BLAT (any version)

###Python requirements

numpy, scipy, biopython, gzip, screed, statsmodel (optional)

## Installation

Clone the git repository

   ```bash
   $> git clone https://github.com/jccastrog/imGLAD
   ```

You can also download the zip file from the GitHub site [https://github.com/jccastrog/imGLAD](https://github.com/jccastrog/imGLAD)

## Building models for target genomes

Once you have installed imGLAD you can use [fitModel](./fitModel.py) to create a model of the target genome you want to detect.

The automatic training generates reads form a randomly selected number of genomes (default is 200 genomes) from RefSeq (Pruitt et al., 2004), and builds in-silico-generated datasets of about 1 million reads each. Simulated reads from the target genome(s) are then generated in a similar way, and added to the former datasets, at different abundances, in order to create the positive datasets. Reads from the target genome(s) are omitted for the construction of negative datasets. All other genomes used to create the datasets are sampled in equal proportions (i.e., even richness).

## Estimating the probability of presence for your samples

Once the logistic model has been built, sequencing breadth can be used to reliably predict the probability of presence of the target genome in any number of query metagenomic datasets, using [probEstimate](./probEstimate). 
