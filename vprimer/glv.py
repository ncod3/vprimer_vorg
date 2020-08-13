# -*- coding: utf-8 -*-

import sys
import os
import errno
import time

# using class
from vprimer.param import Param
from vprimer.conf import Conf
from vprimer.ref_fasta import RefFasta
from vprimer.outlist import OutList

# https://qiita.com/mimitaro/items/3506a444f325c6f980b2
import configparser
# https://docs.python.org/ja/3/library/configparser.html

##########################################
# allele status
AL_HOMO = 'homo'
AL_HETERO = 'hetero'

##########################################
# pick_mode
# all
# INDEL     InDel
# MNVIND     Multiple nucleotide variants & mini indel
# SNP       Single nucleotide polymorphysms
##########################################
MODE_ALL    = "all"
MODE_INDEL  = "indel"    # indel
MODE_SNPMNV = "snpmnv"   # snp mnv mind
MODE_SNP    = "snp"      # snp

##########################################
# var_type (variation_type)
OutOfRange = "oor"
#: Code for indel allele, includes substitutions of unequal length
INDEL = "indel"
#: Code for single nucleotide variant allele
SNP =   "snp"
#: Code for a multi nucleotide variant allele
MNV =   "mnv"
#: mini_indel 1 =< diff_len min_indel_len
MIND =  "mind"

# long_smpl_side
SAME_LENGTH = -1
# caps not digested side
#NOT_DIGESTED = -1

# SKIP
SKIP_DONT_SKIP = -1
SKIP_SAME_HOMO = 1
SKIP_SAME_HETERO = 2
SKIP_None = 3
SKIP_DIFF_INGROUP = 10

# formtxt COMMENT
COMMENT_nop = '-'
COMMENT_dup = 'dup'
COMMENT_blast = 'blast'
COMMENT_AABC = 'AA_BC'
COMMENT_ABCD = 'AB_CD'

###########################################
# diff pattern
#   [A, B, C, D]
#    | AA     BB     CC     DD     AB     AC     AD     BC     BD     CD
# ---+---------------------------------------------------------------------
# AA |        AA,BB  AA,CC  AA,DD  AA,AB  AA,AC  AA,AD  AA,BC  AA,BD  AA,CD
#    |        hoho   hoho   hoho   hohe_s hohe_s hohe_s hohe_n hohe_n hohe_n
#
# BB | BB,AA         BB,CC  BB,DD  BB,AB  BB,AC  BB,AD  BB,BC  BB,BD  BB,CD
#    | hoho          hoho   hoho   hohe_s hohe_n hohe_n hohe_s hohe_s hohe_n
#
# CC | CC,AA  CC,BB         CC,DD  CC,AB  CC,AC  CC,AD  CC,BC  CC,BD  CC,CD
#    | hoho   hoho          hoho   hohe_n hohe_s hohe_n hohe_s hohe_n hohe_s
#
# DD | DD,AA  DD,BB  DD,CC         DD,AB  DD,AC  DD,AD  DD,BC  DD,BD  DD,CD
#    | hoho   hoho   hoho          hohe_n hohe_n hohe_s hohe_n hohe_s hohe_s
#
# AB | AB,AA  AB,BB  AB,CC  AB,DD         AB,AC  AB,AD  AB,BC  AB,BD  AB,CD
#    | hohe_s hohe_s hohe_n hohe_n        hehe_s hehe_s hehe_s hehe_s hehe_n
#
# AC | AC,AA  AC,BB  AC,CC  AC,DD  AC,AB         AC,AD  AC,BC  AC,BD  AC,CD
#    | hohe_s hohe_n hohe_s hohe_n hehe_s        hehe_s hehe_s hehe_n hehe_s
#
# AD | AD,AA  AD,BB  AD,CC  AD,DD  AD,AB  AD,AC         AD,BC  AD,BD  AD,CD
#    | hohe_s hohe_n hohe_n hohe_s hehe_s hehe_s        hehe_n hehe_s hehe_s
#
# BC | BC,AA  BC,BB  BC,CC  BC,DD  BC,AB  BC,AC  BC,AD         BC,BD  BC,CD
#    | hohe_n hohe_s hohe_s hohe_n hehe_s hehe_s hehe_n        hehe_s hehe_s
#
# BD | BD,AA  BD,BB  BD,CC  BD,DD  BD,AB  BD,AC  BD,AD  BD,BC         BD,CD
#    | hohe_n hohe_s hohe_n hohe_s hehe_s hehe_n hehe_s hehe_s        hehe_s
#
# CD | CD,AA  CD,BB  CD,CC  CD,DD  CD,AB  CD,AC  CD,AD  CD,BC  CD,BD
#    | hohe_n hohe_n hohe_s hohe_s hehe_n hehe_s hehe_s hehe_s hehe_s

# segregation_pattern

segr_ptn_NOP                        = 'nop'
# ./.
segr_ptn_NOT_EXIST_ALLELE           = 'not_exist_allele'

# AA,AA BB,BB CC,CC DD,DD
segr_ptn_SAME_HOMO                  = 'same_homo'

# AB,AB AC,AC AD,AD BC,BC BD,BD CD,CD
segr_ptn_SAME_HETERO                = 'same_hetero'

# AA,BB AA,CC AA,DD BB,CC BB,DD CC,DD
# BB,AA CC,AA DD,AA CC,BB DD,BB DD,CC
segr_ptn_HOMO_HOMO                  = 'hoho_1'

# AA,AB AA,AC AA,AD BB,AB BB,BC BB,BD CC,AC CC,BC CC,CD DD,AD
# DD,BD DD,CD
# AB,AA AC,AA AD,AA AB,BB BC,BB BD,BB AC,CC BC,CC CD,CC AD,DD
# BD,DD CD,DD
segr_ptn_HOMO_HETERO_SHARE          = 'hohe_s1'

# AA,BC AA,BD AA,CD BB,AC BB,AD BB,CD CC,AB CC,AD CC,BD DD,AB
# DD,AC DD,BC
# BC,AA BD,AA CD,AA AC,BB AD,BB CD,BB AB,CC AD,CC BD,CC AB,DD
# AC,DD BC,DD
segr_ptn_HOMO_HETERO_NOT_SHARE      = 'hohe_n2'

# AB,AC AB,AD AB,BC AB,BD AC,AD AC,BC AC,CD AD,BD AD,CD BC,BD
# BC,CD BD,CD
segr_ptn_HETERO_HETERO_SHARE        = 'hehe_s3'

# AB,CD AC,BD AD,BC
segr_ptn_HETERO_HETERO_NOT_SHARE    = 'hehe_n4'

# 使用する制限酵素の最長のrecognition site length
AROUND_SEQ_LEN = 20

def init(prog_name):

    global program_name
    program_name = prog_name

    global param
    param = Param()

    global conf
    conf = Conf()

    global ref
    ref = RefFasta()

    global outlist
    outlist = OutList()

    # read parameter into global environment
    param = param.get_args()

    # read ini file into global environment
    conf = conf.read_ini()
    conf = conf.merge_conf()

    # open_log in glv
    param.open_log()
    ref.open_log()
    outlist.open_log()

