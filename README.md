# imGLAD

imGLAD is a computational tool for detection of bacterial genomes in metagenomic datasets.

## Overview

The software consists of two parts the first part creates a series of metagenomic datasets, the datasets are created in such a way that the target organism is present in half of them and absent in the other half.

##System requirements

Python 2.7 or higher.

ART 2.5.8 or higher.

Either BLAST 2.2.28 (or higher) or BLAT (any version)

###Python requirements

numpy, scipy, biopython, gzip, screed

## Installation

Clone the git repository

   ```bash
   $> git clone https://github.com/jccastrog/imGLAD
   ```

You can also download the zip file from the GitHub site https://github.com/jccastrog/imGLAD(https://github.com/jccastrog/imGLAD)

## Building models for target genomes

Once you have installed imGLAD you can use fitModel.py(./fitModel.py) to create a model of the target genome you want to detect. 
