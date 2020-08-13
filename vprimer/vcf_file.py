#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import errno
import time
import re

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
log = LogConf.open_log(__name__)

import vcfpy
import subprocess as sbp


class VCFFile(object):
    """
    """

    def __init__(self):
        pass

    #-------------------------------------------
    @classmethod
    def prepare_vcf(cls):

        # user's vcf: convert relative path to absolute path based on cwd
        if glv.conf.vcf.startswith('/'):
            # originally absolute path
            glv.conf.vcf_file_user = glv.conf.vcf
        else:
            # cwd + relative path
            glv.conf.vcf_file_user = "{}/{}".format(
                glv.conf.cwd,
                glv.conf.vcf)

        log.info("glv.conf.vcf_file_user {}".format(glv.conf.vcf_file_user))

        # vcf_file_user: user asigned vcf in ini file
        if os.path.isfile(glv.conf.vcf_file_user):
            log.info("{} found.".format(glv.conf.vcf_file_user))
        else:
            log.info("{} not found. exit.".format(glv.conf.vcf_file_user))
            sys.exit(1)

        # vcf_file as slink org_slink.gz
        basename_user = os.path.basename(glv.conf.vcf_file_user)
        glv.conf.vcf_file_slink_system = "{}/{}.{}".format(
            glv.conf.ref_dir, basename_user, 'org_slink.gz')

        #log.debug("glv.conf.ref_dir={}".format(glv.conf.ref_dir))
        #log.debug("glv.conf.log_dir={}".format(glv.conf.log_dir))
        #log.debug("glv.conf.out_bak_dir={}".format(glv.conf.out_bak_dir))
        #log.debug("glv.conf.log_dir={}".format(glv.conf.log_dir))

        # symlink glv.conf.vcf_file_user to glv.conf.vcf_file
        if os.path.isfile(glv.conf.vcf_file_slink_system):
            log.info("{} exist.".format(glv.conf.vcf_file_slink_system))
        else:
            os.symlink(glv.conf.vcf_file_user, glv.conf.vcf_file_slink_system)
            log.info("os.symlink {} {}.".format(
                glv.conf.vcf_file_user, glv.conf.vcf_file_slink_system))

        # gtonly.gz
        glv.conf.vcf_file = re.sub(
            r"\.org_slink.gz$", "GTonly.gz", glv.conf.vcf_file_slink_system)

        # glv.conf.vcf_file
        if os.path.isfile(glv.conf.vcf_file):
            log.info("{} exist.".format(glv.conf.vcf_file))
        else:
            cmd1 = "{} annotate --threads {} -O z -x ^FMT/GT -o {} {}".format(
                'bcftools',
                glv.conf.parallele_full_thread,
                glv.conf.vcf_file,
                glv.conf.vcf_file_slink_system)


        # if bcftools 1.10.2, it have --force option.
        # if 1.9.0, if we have an error
        # No matching tag in -x ^FMT/GT
        # we will not worry, continue

            err_str = utl.try_exec_error(cmd1)
            if err_str == '' or \
                err_str.startswith('No matching tag in'):
                log.error("we will go ahead if bcftools says, {}.".format(
                    err_str))

                os.remove(glv.conf.vcf_file)

                os.symlink(glv.conf.vcf_file_slink_system, glv.conf.vcf_file)
                log.info("os.symlink {} {}.".format(
                    glv.conf.vcf_file_slink_system, glv.conf.vcf_file))

            else:
                log.error("{}.".format(err_str))
                sys.exit(1)

            # make tbi
            utl.tabix(glv.conf.vcf_file)


        # sample names to file
        sample_name_file = "{}.sample_name.txt".format(glv.conf.vcf_file)

        if os.path.isfile(sample_name_file):
            cls._read_sample_name_list(sample_name_file)
            log.info("exist sample_name_file {}".format(sample_name_file))

        else:
            sample_names = cls._get_sample_name_list()
            log.info("not exist {}".format(sample_name_file))

            with open(sample_name_file, mode='w') as f:
                # write list
                f.write("{}\n".format("\n".join(sample_names)))
            cls._read_sample_name_list(sample_name_file)


    @classmethod
    def _get_sample_name_list(cls):

        basename_fullname = list()
        reader = vcfpy.Reader.from_path(glv.conf.vcf_file)
        # list to text

        basename_user = os.path.basename(glv.conf.vcf_file_user)

        for (sno, sample_fullname) in enumerate(
            reader.header.samples.names, 1):

            sample_basename = os.path.basename(sample_fullname)
            basename_fullname.append("{}\t{}\t{}".format(
                sno, sample_basename, sample_fullname))

            log.debug("{} {}".format(sample_basename, sample_fullname))

        #log.debug("sample_names {}".format(sample_names))
        
        # no  basename            fullname
        # 1   sample1_sorted.bam  /home/user/sample1_sorted.bam
        # 2   sample2_sorted.bam  /home/user/sample2_sorted.bam

        #glv.conf.vcf_sample_name[sample_basename] = sample_fullname

        # glv.conf.vcf_sample_name = {
        #   'sample1_sorted.bam' : '/home/user/sample1_sorted.bam',
        #   'sample2_sorted.bam' : '/home/user/sample2_sorted.bam' }

        return basename_fullname


    @classmethod
    def _read_sample_name_list(cls, sample_name_file):

        with open(sample_name_file, mode='r') as r:
            for liner in r:
                r_line = liner.strip()

                sno, sample_basename, sample_fullname = \
                    r_line.split('\t')
                cls._check_sample_name(sample_basename, sample_fullname)


    @classmethod
    def _check_sample_name(cls, sample_basename, sample_fullname):

        #log.debug("basename={}, fullname={}".format(
        #    sample_basename, sample_fullname))

        glv.conf.vcf_basename_to_fullname[sample_basename] = \
            sample_fullname

        # 1) is basename have nickname?
        if sample_basename in glv.conf.basename_to_nickname:
            
            nick_name = glv.conf.basename_to_nickname[sample_basename]
            #log.debug("basename={}, nick_name={}".format(
            #    sample_basename, nick_name))

            #if nick_name in glv.conf.nickname_to_basename[nick_name]:
            if nick_name in glv.conf.nickname_to_basename:
                glv.conf.vcf_nickname_to_fullname[nick_name] = \
                    sample_fullname
                #log.debug("nick_name={}, fullname={}".format(
                #    nick_name, sample_fullname))

            else:
                log.error("nick_name {} has not in vcf.".format(
                    nick_name))
                sys.exit(1)


