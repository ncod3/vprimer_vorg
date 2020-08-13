#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import glob
import pandas as pd
import time
import datetime

import logging
import logging.config

# global configuration
import vprimer.glv as glv

import subprocess as sbp
from subprocess import PIPE


def start_log():
    ''' from conf.py _set_paths
    '''

    logging.config.dictConfig(glv.conf.log.config)
    global log
    log = logging.getLogger(__name__)
    log.info("logging start {}".format(__name__))


def save_to_tmpfile(file_path, can_log = True):
    """
    """

    ret = False
    if os.path.isfile(file_path):
        # /a/b/c.txt
        # /a/b/bak/c.txt
        dirname_file = os.path.dirname(file_path)
        basename_file = os.path.basename(file_path)

        file_bak_path = "{}/{}".format(
            glv.conf.out_bak_dir,
            basename_file)

        ret = True
        ts = time.time()
        new_file_path = "{}.{}.bak".format(file_bak_path, ts)

        os.rename(file_path, new_file_path)

        if can_log:
            log.info("{} exist. mv to {}".format(
                file_path, new_file_path))

    return ret


def write_df_to_csv(dataframe, out_file, index=True, force=True):
    """
    """

    file_name = out_file
    moved = False
    if force == True:
        moved = save_to_tmpfile(out_file)

    dataframe.to_csv(out_file, sep='\t', index=index)


def progress_check(now_progress):

    stat = False    # False if don't do this progress
    param_progress = glv.conf.progress

    log.info("now_progress={} param_progress={}".format(
        now_progress, param_progress))

    #log.debug("now_progress={} {} param_progress={} {}".format(
    #    now_progress,
    #    now_progress_no,
    #    param_progress,
    #    param_progress_no))

    if param_progress == 'all':
        stat = True

    else:
        now_progress_no = glv.outlist.outf_prefix[now_progress]['no']
        param_progress_no = glv.outlist.outf_prefix[param_progress]['no']
        if now_progress_no >= param_progress_no:
            stat = True

    return stat

def stop(now_progress):

    if glv.conf.stop == 'no':
        return

    now_progress_no = glv.outlist.outf_prefix[now_progress]['no']
    param_stop_no = glv.outlist.outf_prefix[glv.conf.stop]['no']
    if now_progress_no >= param_stop_no:
        log.info("stop {}".format(glv.conf.stop))
        sys.exit(1)


def is_my_pick_mode(var_type, pick_mode):

    if var_type == glv.OutOfRange:
        return False

    elif pick_mode == glv.MODE_ALL or \
       (pick_mode == glv.MODE_INDEL  and  var_type == glv.INDEL) or \
       (pick_mode == glv.MODE_SNP    and  var_type == glv.SNP) or \
       (pick_mode != glv.MODE_INDEL  and  var_type != glv.INDEL):
        return True

    else:
        return False


def check_for_files(filepath):

    # filepath is pattern
    fobj_list = list()

    for filepath_object in glob.glob(filepath):

        if os.path.isfile(filepath_object):
            fobj_list.append(filepath_object)

    return sorted(fobj_list)


def is_None(s0_0, s0_1, s1_0, s1_1):
    return s0_0 == None or s0_1 == None or \
           s1_0 == None or s1_1 == None


def is_homo_homo(s0_0, s0_1, s1_0, s1_1):
           #  1 == 2           3 == 4
    return s0_0 == s0_1 and s1_0 == s1_1


def is_homo_hetero(s0_0, s0_1, s1_0, s1_1):
           #  1 == 2           3 != 4
    return s0_0 == s0_1 and s1_0 != s1_1 or \
           s0_0 != s0_1 and s1_0 == s1_1
           #  1 != 2           3 == 4


def is_hetero_hetero(s0_0, s0_1, s1_0, s1_1):
           #  1 != 2           3 != 4
    return s0_0 != s0_1 and s1_0 != s1_1


def is_same_homo(s0_0, s0_1, s1_0, s1_1):
    return is_homo_homo(s0_0, s0_1, s1_0, s1_1) and \
           s0_0 == s1_0
           #  1 == 3


def is_same_hetero(s0_0, s0_1, s1_0, s1_1):
    return is_hetero_hetero(s0_0, s0_1, s1_0, s1_1) and \
           s0_0 == s1_0 and s0_1 == s1_1
           #  1 == 3           2 == 4

def is_share(s0_0, s0_1, s1_0, s1_1):
           #  1 == 3          1 == 4
    return s0_0 == s1_0 or s0_0 == s1_1 or \
           s0_1 == s1_0 or s0_1 == s1_1
           #  2 == 3          2 == 4

def strip_hash_comment(line):
    return line.split('#')[0].strip()


def elapsed_time(now_time, start_time):

    elapsed_time = now_time - start_time
    #3:46:11.931354
    return "elapsed_time {}".format(
        datetime.timedelta(seconds=elapsed_time))


def try_exec_error(cmd):
    # https://qiita.com/HidKamiya/items/e192a55371a2961ca8a4

    err_str = ''

    try:
        log.info("do {}".format(cmd))

        sbp.run(cmd,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
            shell=True,
            check=True)

    except sbp.CalledProcessError as e:
        err_str = e.stderr

    return err_str


def try_exec(cmd):

    try:
        log.info("do {}".format(cmd))

        sbp.run(cmd,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
            shell=True,
            check=True)

    except sbp.CalledProcessError as e:
        log.error("{}.".format(e.stderr))
        sys.exit(1)


def tabix(vcf_file):

    #-f -p vcf
    cmd1 = "{} {} {}".format(
        'tabix',
        '-f -p vcf',
        vcf_file)

    try_exec(cmd1)


def sort_file(
        proc, distin_dict, out_file_name,
        nm_chrom, nm_pos, nm_order, n):

    # sort command option
    if n == 'number':
        n = 'n'
    else:
        n = ''


    # cmd.sort
    hdr_dict = distin_dict[proc]['hdr_dict']
    out_txt_file = distin_dict[proc]['out_path']
    sorted_file = "{}.sorted".format(out_txt_file)

    col_chrom = hdr_dict[nm_chrom]
    col_pos = hdr_dict[nm_pos]
    col_order = hdr_dict[nm_order]

    cmd_sort = "{} -k {},{} -k {},{}n -k {},{}{} {} > {}".format(
        'sort',
        col_chrom, col_chrom,
        col_pos, col_pos,
        col_order, col_order,
        n,
        out_txt_file,
        sorted_file)

    try_exec(cmd_sort)

    # rm file
    log.info("remove {}".format(out_txt_file))
    os.remove(out_txt_file)


    # make header.txt
    out_txt_header = "{}.header_txt".format(out_txt_file)
    with open(out_txt_header, mode='w') as f:
        f.write("{}\n".format(distin_dict[proc]['hdr_text']))

    # cat header file
    cmd_cat = "{} {} {} > {}".format(
        'cat',
        out_txt_header,
        sorted_file,
        out_txt_file)

    try_exec(cmd_cat)

    # rm header sorted
    os.remove(out_txt_header)
    os.remove(sorted_file)




