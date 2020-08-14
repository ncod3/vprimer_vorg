#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# https://genome.sph.umich.edu/wiki/Variant_classification

import sys
import os
import errno
import time

#import logging
#import vprimer.logging_config
#log = logging.getLogger(__name__)

# global variables
import vprimer.glv as glv
import vprimer.utils as utl
from vprimer.logging_config import LogConf

#--- read dir and log file set
glv.init(sys.argv[0])
log = LogConf.open_log(__name__)
log.info("logging start {}".format(__name__))

# using class
from vprimer.vcf_file import VCFFile
from vprimer.variant import Variant
from vprimer.marker import Marker
from vprimer.primer import Primer
from vprimer.enzyme import Enzyme
from vprimer.formtxt import FormTxt


def main():

    start = time.time()
    log.info('program started')

    # run
    vpr = VPrimer()
    vpr.run()

    log.info("program finished {}".format(
        utl.elapsed_time(time.time(), start)))


class VPrimer(object):

    def __init__(self):

        self.vcf = VCFFile()
        self.variant = Variant()
        self.marker = Marker()
        self.primer = Primer()
        self.formtxt = FormTxt()

    def run(self):

        # prepare
        self.prepare()
        utl.stop('prepare')

        # pick_variant
        self.variant.pick_variant()
        utl.stop('variant')

        # design_marker 
        self.marker.design_marker()
        utl.stop('marker')

        # primer3 input file
        self.primer.construct_primer()
        utl.stop('primer')

        # output as txt
        self.formtxt.format_text()


    def prepare(self):

        # read reference into global environment
        glv.ref = glv.ref.prepare_ref()

        # only prepare vcf file @classmethod
        VCFFile.prepare_vcf()

        # prepare ensyme files @classmethod
        Enzyme.prepare_enzyme()

        # prepare output text by glv.conf.distin_g_list
        glv.outlist.prepare_distin_files()


if __name__ == '__main__':
    main()

