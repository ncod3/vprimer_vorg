import sys
import os
import errno
import time

import logging
log = logging.getLogger(__name__)

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

class AlleleSelect(object):

    def __init__(self):

        self.lines = list()
        self.var_types = list()

        self.top_smpl_list = list()
        self.gr_list = list()

        self.chrom = ''
        self.pos = 0
        self.alt_vc = 0
        self.vseq_ano = list()
        self.vseq_ano_str = ''
        self.len_ano = list()

        # indel size threshold
        self.min_indel_len = glv.conf.min_indel_len
        self.max_indel_len = glv.conf.max_indel_len

        # ano_Group_alleleIdx
        self.g0x0 = -1
        self.g0x1 = -1
        self.g1x0 = -1
        self.g1x1 = -1

        # some list
        self.ano__gno_aix = list()  # from record

        self.rgt__gno_aix = list()
        self.gti__gno_aix = list()
        self.ali__gno_aix = list()
        self.len__gno_aix = list()

        self.altype_gno = list()
        self.alint_gno = list()

        self.segr_ptn = glv.segr_ptn_NOP
        self.diff_allele_set = list()
        self.diff_allele_cnt = 0


    def select_diff_allele(self, record, top_smpl_list, gr_list):

        self.top_smpl_list = top_smpl_list  # ['Akenohoshi', 'Akitakomachi']
        self.gr_list = gr_list              # ['gAkenohoshi', 'gAkitakomachi']

        # get basic record info & ano__gno_aix
        self._get_record_info(record)

        # get basic sample's allele combination info
        self._get_allele_info()

        # get allele combination for marker among four alleles in two sample
        self._get_segregation_pattern()

        # for print
        for no, diff_allele in enumerate(self.diff_allele_set, 1):

            var_type, line = self._construct_var_line(no, diff_allele)
            self.lines.append(line)
            self.var_types.append(var_type)


    def _construct_var_line(self, no, diff_allele):

        line_list = list()

        # g0
        g0_ano = diff_allele[0]
        g0_len = self.len_ano[g0_ano]
        g0_name = self.gr_list[0]

        # g1
        g1_ano = diff_allele[1]
        g1_len = self.len_ano[g1_ano]
        g1_name = self.gr_list[1]

        # gts_segr_lens
        gts_str, lens_str = self._get_gts_lens_str()
        gts_segr_lens = "{},{},{}".format(
            gts_str, self.segr_ptn, lens_str)

        # var_type, longest_gno, longest_len, diff_len
        var_type, longest_gno, longest_len, diff_len = \
            self._get_variant_type(g0_ano, g1_ano)

        targ_grp = "{},{}".format(g0_name, g1_name)
        targ_ano = "{},{}".format(g0_ano, g1_ano)
        set_n = "{}/{}".format(no, self.diff_allele_cnt)
        len_g0g1_dif_long = "{},{},{},{}".format(
            g0_len, g1_len, diff_len, longest_gno)

        # ---------------------------------------
        line_list += [self.chrom]
        line_list += [self.pos]
        line_list += [targ_grp]
        line_list += [gts_segr_lens]
        line_list += [set_n]
        line_list += [targ_ano]
        line_list += [len_g0g1_dif_long]
        line_list += [var_type]
        line_list += [self.vseq_ano_str.upper()]

        return var_type, '\t'.join(map(str, line_list))


    def _get_gts_lens_str(self):

        # gt_segr_len  genotypes and segregation type and allele len
        gtsegr_gtlist = list()
        gtsegr_lenlist = list()

        for gno in range(2):
            # for gt_segr_len
            gtsegr_gtlist.append("".join(map(str, self.gti__gno_aix[gno])))
            gtsegr_lenlist.append(".".join(map(str, self.len__gno_aix[gno])))

        gts_str = "/".join(map(str, gtsegr_gtlist))
        lens_str = "/".join(map(str, gtsegr_lenlist))

        return gts_str, lens_str


    def _get_variant_type(self, g0_ano, g1_ano):

        # glv.SNP
        # glv.MNV
        # glv.MIND
        # glv.INDEL
        # glv.OutOfRange

        var_type = glv.OutOfRange 

        g0_len = self.len_ano[g0_ano]
        g1_len = self.len_ano[g1_ano]

        longest_gno = glv.SAME_LENGTH
        longest_len = g0_len

        if g0_len < g1_len:
            longest_gno = 1
            longest_len = g1_len
        elif g0_len > g1_len:
            longest_gno = 0
            longest_len = g0_len

        diff_len = abs(g0_len - g1_len)

        if diff_len == 0 and g0_len == 1:
            var_type = glv.SNP

        elif diff_len == 0 and longest_len > 1:
            self.var_type = glv.MNV

        elif diff_len < self.min_indel_len:
            var_type = glv.MIND

        elif self.min_indel_len <= diff_len and \
            diff_len <= self.max_indel_len:
            var_type = glv.INDEL

        else:
            var_type = glv.OutOfRange

        return var_type, longest_gno, longest_len, diff_len
        

    def _get_segregation_pattern(self):

        # -- already removed
        # segr_ptn_NOT_EXIST_ALLELE           = 'not_exist_allele'
        # segr_ptn_SAME_HOMO                  = 'same_homo'
        # segr_ptn_SAME_HETERO                = 'same_hetero'

        # -- check now
        # segr_ptn_HOMO_HOMO                  = 'hoho'
        # segr_ptn_HOMO_HETERO_SHARE          = 'hohe_s'
        # segr_ptn_HOMO_HETERO_NOT_SHARE      = 'hohe_n'
        # segr_ptn_HETERO_HETERO_SHARE        = 'hehe_s'
        # segr_ptn_HETERO_HETERO_NOT_SHARE    = 'hehe_n'

        # self.g0x0, self.g0x1, self.g1x0, self.g1x1
        #log.debug("g0=[{}, {}], g1=[{}, {}]".format(
        #    self.g0x0, self.g0x1, self.g1x0, self.g1x1))

        #             set_n
        # 1.homo vs homo
        #   hoho        1     00/11 0,1
        # 2.homo vs hetero
        #   hohe_s      1     00/01 0,1
        #   hohe_n      2     00/12 0,1 0,2
        # 3.hetero vs hetero
        #   hehe_s      3     01/02 0,2 1,0 1,2
        #   hehe_n      4     01/23 0,2 0,3 1,2 1,3

        # 1.homo vs homo
        if utl.is_homo_homo(self.g0x0, self.g0x1, self.g1x0, self.g1x1):
            # AA,BB
            #   hoho        1     00/11 0,1
            self.segr_ptn = glv.segr_ptn_HOMO_HOMO
            #                           [      AA0,       BB0]
            self.diff_allele_set.append([self.g0x0, self.g1x0])

        # 2.homo vs hetero
        elif utl.is_homo_hetero(self.g0x0, self.g0x1, self.g1x0, self.g1x1):

            if utl.is_share(self.g0x0, self.g0x1, self.g1x0, self.g1x1):
                # AA,AB
                #   hohe_s      1     00/01 0,1
                self.segr_ptn = glv.segr_ptn_HOMO_HETERO_SHARE

                if self.altype_gno[0] == glv.AL_HETERO:
                    # AB,AA(BB)
                    if self.g0x0 != self.g1x0:
                        # AB,BB                     [      AB0,       BB0]
                        self.diff_allele_set.append([self.g0x0, self.g1x0])
                    else:
                        # AB,AA                     [      AB1,       AA0]
                        self.diff_allele_set.append([self.g0x1, self.g1x0])

                else:
                    # AA(BB),AB
                    if self.g0x0 != self.g1x0:
                        # BB,AB                     [      BB0,       AB0]
                        self.diff_allele_set.append([self.g0x0, self.g1x0])
                    else:
                        # AA,AB                     [      AA0,       AB1]
                        self.diff_allele_set.append([self.g0x0, self.g1x1])

            else:
                # AA,BC
                #   hohe_n      2     00/12 0,1 0,2
                self.segr_ptn = glv.segr_ptn_HOMO_HETERO_NOT_SHARE

                if self.altype_gno[0] == glv.AL_HETERO:
                    # BC,AA -> [B,A] [C,A]
                    self.diff_allele_set.append([self.g0x0, self.g1x0])
                    self.diff_allele_set.append([self.g0x1, self.g1x0])

                else:
                    # AA,BC -> [A,B] [A,C]
                    self.diff_allele_set.append([self.g0x0, self.g1x0])
                    self.diff_allele_set.append([self.g0x0, self.g1x1])

        # 3.hetero vs hetero
        else:

            if utl.is_share(self.g0x0, self.g0x1, self.g1x0, self.g1x1):
                # AB,AC
                #   hehe_s      3     01/02 0,2 1,0 1,2
                self.segr_ptn = glv.segr_ptn_HETERO_HETERO_SHARE


                # 01,02
                if self.g0x0 == self.g1x0:
                    # 0/2, 1/0, 1/2
                    #self.diff_allele_set.append([self.g0x0, self.g1x0])
                    self.diff_allele_set.append([self.g0x0, self.g1x1])
                    self.diff_allele_set.append([self.g0x1, self.g1x0])
                    self.diff_allele_set.append([self.g0x1, self.g1x1])

                # 01,20
                elif self.g0x0 == self.g1x1:
                    # 0/2, 1,2, 1,0
                    self.diff_allele_set.append([self.g0x0, self.g1x0])
                    #self.diff_allele_set.append([self.g0x0, self.g1x1])
                    self.diff_allele_set.append([self.g0x1, self.g1x0])
                    self.diff_allele_set.append([self.g0x1, self.g1x1])

                # 10,02
                elif self.g0x1 == self.g1x0:
                    # 1,0, 1,2, 0,2
                    self.diff_allele_set.append([self.g0x0, self.g1x0])
                    self.diff_allele_set.append([self.g0x0, self.g1x1])
                    #self.diff_allele_set.append([self.g0x1, self.g1x0])
                    self.diff_allele_set.append([self.g0x1, self.g1x1])

                # 10,20
                elif self.g0x1 == self.g1x1:
                    # 1,2, 1,0, 0,2
                    self.diff_allele_set.append([self.g0x0, self.g1x0])
                    self.diff_allele_set.append([self.g0x0, self.g1x1])
                    self.diff_allele_set.append([self.g0x1, self.g1x0])
                    #self.diff_allele_set.append([self.g0x1, self.g1x1])

            else:
                # AB,CD
                #   hehe_n      4     01/23 0,2 0,3 1,2 1,3
                self.segr_ptn = glv.segr_ptn_HETERO_HETERO_NOT_SHARE
                # all combination 4 type
                #                                   A.         C.
                self.diff_allele_set.append([self.g0x0, self.g1x0])
                #                                   A.         .D
                self.diff_allele_set.append([self.g0x0, self.g1x1])
                #                                   .B         C.
                self.diff_allele_set.append([self.g0x1, self.g1x0])
                #                                   .B         .D
                self.diff_allele_set.append([self.g0x1, self.g1x1])

        self.diff_allele_cnt = len(self.diff_allele_set)


    def _get_allele_info(self):

        # gno 0..1 aix 0..1
        for gno in range(2):

            rgt_aix = list()    # raw_gt
            gti_aix = list()    # None->-1
            ali_aix = list()    # +1
            len_aix = list()

            for aix in range(2):
                raw_gt = self.ano__gno_aix[gno][aix]
                if raw_gt == None:  # ./.
                    rgt_aix.append(',')
                    ali = -1 + 1
                    lena = 0
                    gti_aix.append(-1)
                else:
                    ali = raw_gt + 1
                    rgt_aix.append(raw_gt)
                    gti_aix.append(raw_gt)
                    lena = self.len_ano[raw_gt]

                ali_aix.append(ali)
                len_aix.append(lena)

            # row genotype [gno][aix]
            self.rgt__gno_aix.append(rgt_aix)

            # genotype integer [gno][aix] (None == 0)
            self.gti__gno_aix.append(gti_aix)

            # allele integer [gno][aix]
            self.ali__gno_aix.append(ali_aix)

            # allele len [gno][aix]
            self.len__gno_aix.append(len_aix)

            # self.altype_gno
            if self.ano__gno_aix[gno][0] == self.ano__gno_aix[gno][1]:
                self.altype_gno.append(glv.AL_HOMO)
            else:
                self.altype_gno.append(glv.AL_HETERO)

            # allele int
            self.alint_gno.append(
                (self.ali__gno_aix[gno][0] * 10) + self.ali__gno_aix[gno][1])


    def _get_record_info(self, record):

        self.chrom = record.CHROM
        self.pos = record.POS
        self.alt_vc = len(record.ALT)   # alt variant count
        # get vseq from ano
        self.vseq_ano = [record.REF]
        [self.vseq_ano.append(alt.value) for alt in record.ALT]
        # str
        self.vseq_ano_str = ','.join(map(str, self.vseq_ano))
        # len
        [self.len_ano.append(len(sqa)) for sqa in self.vseq_ano]
        # get ano from gno x aix
        sample0 = self.top_smpl_list[0]
        sample1 = self.top_smpl_list[1]

        # get ano from record by sample name
        self.g0x0, self.g0x1, self.g1x0, self.g1x1 = \
            AlleleSelect.record_call_for_sample(record, sample0, sample1)

        # ano__gno_aix
        self.ano__gno_aix = [[self.g0x0, self.g0x1], [self.g1x0, self.g1x1]]

    @classmethod
    def record_call_for_sample(cls, record, sample0, sample1):

        fullname0 = cls._get_sample_fullname(sample0)
        fullname1 = cls._get_sample_fullname(sample1)

        # for REF 20200708
        if sample0 == 'REF':
            s0_0 = 0
            s0_1 = 0
            s1_0 = record.call_for_sample[fullname1].gt_alleles[0]
            s1_1 = record.call_for_sample[fullname1].gt_alleles[1]

        elif sample1 == 'REF':
            s0_0 = record.call_for_sample[fullname0].gt_alleles[0]
            s0_1 = record.call_for_sample[fullname0].gt_alleles[1]
            s1_0 = 0
            s1_1 = 0

        else:
            s0_0 = record.call_for_sample[fullname0].gt_alleles[0]
            s0_1 = record.call_for_sample[fullname0].gt_alleles[1]
            s1_0 = record.call_for_sample[fullname1].gt_alleles[0]
            s1_1 = record.call_for_sample[fullname1].gt_alleles[1]

        # int
        return s0_0, s0_1, s1_0, s1_1


    @classmethod
    def _get_sample_fullname(cls, name):

        fullname = ''

        #log.debug("n2f={}".format(glv.conf.vcf_nickname_to_fullname))
        #log.debug("b2f={}".format(glv.conf.vcf_basename_to_fullname))

        if name in glv.conf.vcf_nickname_to_fullname:
            fullname = glv.conf.vcf_nickname_to_fullname[name]
            #log.debug("1 name={} fullname={}".format(name, fullname))

        elif name in glv.conf.vcf_basename_to_fullname:
            fullname = glv.conf.vcf_basename_to_fullname[name]
            #log.debug("2 name={} fullname={}".format(name, fullname))

        else:
            fullname = name
            #log.debug("3 name={} fullname={}".format(name, fullname))

        return fullname
    
        

