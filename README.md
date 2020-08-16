# Vprimer

vprimer

## Features

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
blastn
~~~

## Installation
~~~
if you have not yet made vprimer environment on conda, you can make it.
$ conda create -n run_vprimer python=3.7
$ source activate run_vprimer
(run_vprimer) pip install git+https://github.com/ncod3/vprimer


If you have not yet made vprimer environment on conda, you can make it.
$ conda create -n run_vprimer python=3.7
$ source activate run_vprimer
$ pip install git+https://github.com/ncod3/vprimer

If you want to uninstall vprimer,
$ pip uninstall vprimer

~~~

## Getting Started

~~~
If you use minimum cpu 2, it will take only about 5 minutes.
$ vprimer -c ini_yam_exam15.ini -t 2

If you can use more cpu, for example 10, it will take only about 2 minutes.
$ vprimer -c ini_yam_exam15.ini -t 10
~~~

## Usage

## Note

## Authors
- Satoshi Natsume

See also the list of contributors who participated in this project.

## Licence
This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgements

