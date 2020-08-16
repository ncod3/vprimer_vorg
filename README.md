# Vprimer

vprimer

## Features

Test data can download through git clone
- https://github.com/ncod3/data_vprimer

You can make vcf using svaba from several bam and reference.
- https://github.com/ncod3/indelvcf
 

## Contents

## Requirement

~~~
python 3.7
pandas
vcfpy
pysam
joblib
Bio

bcftools
tabix
primer3
blastn
~~~

## Installation
~~~
if you have not yet made vprimer environment on conda, you can make it.

$ conda create -n run_vprimer python=3.7
$ source activate run_vprimer
$ (run_vprimer) pip install git+https://github.com/ncod3/vprimer

We recommend that you start by uninstalling vprimer as it will be updated frequently.

$ (run_vprimer) pip uninstall vprimer

Even if vprimer is not installed, the following message will be output.
WARNING: Skipping vprimer as it is not installed.

$ (run_vprimer) conda install -c anaconda pandas
$ (run_vprimer) conda install -c bioconda vcfpy
$ (run_vprimer) conda install -c bioconda pysam 
$ (run_vprimer) conda install -c anaconda joblib 
$ (run_vprimer) conda install -c anaconda biopython

$ (run_vprimer) conda install -c bioconda bcftools
$ (run_vprimer) conda install -c bioconda tabix
$ (run_vprimer) conda install -c bioconda primer3
$ (run_vprimer) conda install -c bioconda blast

~~~

## Getting Started

~~~
Get the sample data.

$ cd (your working directory)
$ git clone https://github.com/ncod3/data_vprimer
$ mv data_vprimer/ini* .

If you use minimum cpu 2, it will take only about 5 minutes.
$ vprimer -c ini_vprimer_yam_exam15.ini -t 2

If you can use more cpu, for example 10, it will take only about 2 minutes.
$ vprimer -c ini_vprimer_yam_exam15.ini -t 10
~~~

## Usage

If your ini file has a correct vcf file name, you can confirm the sample name in vcf file.
This command only print information, don't touch any thing else.

~~~
(run_vprimer) vprimer -c ini_vprimer_yam_exam15.ini -m
~~~

## Note

## Authors
- Satoshi Natsume s-natsume@ibrc.or.jp

See also the list of contributors who participated in this project.

## Licence
This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgements

