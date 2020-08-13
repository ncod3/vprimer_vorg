# -*- coding: utf-8 -*-

import sys
import os
import errno
import time

# global configuration
import vprimer.glv as glv

import argparse
from multiprocessing import Pool
import multiprocessing as multi
from vprimer.__init__ import __version__
from vprimer.logging_config import LogConf


class Param(object):

    def __init__(self):
        ''' can't write log directly
        '''

    def open_log(self):

        global log
        log = LogConf.open_log(__name__)


    def get_args(self):

        parser = self._vprimer_options()

        if len(sys.argv) == 1:
            self.p = parser.parse_args(['-h'])
        else:
            self.p = parser.parse_args()

        #log.info("{}".format(param))

        return self

    def _vprimer_options(self):

        # https://docs.python.org/ja/3/library/argparse.html
        parser = argparse.ArgumentParser(
            description='vprimer version {}'.format(__version__),
            formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('vprimer ...\n')

        # set options
        parser.add_argument('-V',
                            '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        parser.add_argument('-c',
                            '--config',
                            action='store',
                            required=True,
                            type=str, metavar='',
                            help=".")

        parser.add_argument('-r',
                            '--ref',
                            action='store', #required=True,
                            type=str, metavar='',
                            help=".")

        parser.add_argument('-v',
                            '--vcf',
                            action='store', #required=True,
                            type=str, metavar='',
                            help=".")

        parser.add_argument('-t',
                            '--thread',
                            action='store', #required=True,
                            type=int, metavar='',
                            help=".")

        parser.add_argument('-o',
                            '--out_dir',
                            action='store', #required=True,
                            type=str, metavar='',
                            help=".")

        parser.add_argument('-j',
                            '--joblib_threading',
                            action='store', #required=True,
                            type=str, metavar='',
                            help=".")

        parser.add_argument('-p',
                            '--progress',
                            action='store', #required=True,
                            type=str, metavar='',
                            help=".")

        parser.add_argument('-s',
                            '--stop',
                            action='store', #required=True,
                            type=str, metavar='',
                            help=".")


        #---------------------------------------
        parser.add_argument('-n',
                            '--min_indel_len',
                            action='store', #required=True,
                            type=int, metavar='',
                            help=".")

        parser.add_argument('-x',
                            '--max_indel_len',
                            action='store', #required=True,
                            type=int, metavar='',
                            help=".")

        parser.add_argument('-N',
                            '--min_product_size',
                            action='store', #required=True,
                            type=int, metavar='',
                            help=".")

        parser.add_argument('-X',
                            '--max_product_size',
                            action='store', #required=True,
                            type=int, metavar='',
                            help=".")

        return parser

