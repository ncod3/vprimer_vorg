# -*- coding: utf-8 *-*

import sys
import os
import errno
import time

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
log = LogConf.open_log(__name__)

import vcfpy

from vprimer.allele_select import AlleleSelect

class Variant(object):

    def __init__(self):
        pass


    def pick_variant(self):
        """
        """

        # open vcf through vcfpy
        reader = vcfpy.Reader.from_path(glv.conf.vcf_file)

        # for each distinguish_groups
        for distin_dict in glv.outlist.distin_files:

            #log.debug("distin_dict={}.".format(distin_dict))

            reg = distin_dict['region']
            if reg == '':
                reg = 'whole'
                vcf_ittr = reader
            else:
                vcf_ittr = reader.fetch(glv.conf.regions_dict[reg]['reg'])

            # progress check
            if utl.progress_check('variant') == False:
                log.info("progress={} so skip variant.".format(
                    glv.conf.progress))
                continue

            else:
                log.info("Start processing {}".format('variant'))
                # iterate vcf, using samtools library vcfpy
                self._iterate_vcf(vcf_ittr, distin_dict, reg)


    def _iterate_vcf(self, vcf_ittr, distin_dict, reg):
        """
        """

        pick_mode = distin_dict['pick_mode']
        # 辞書のキーが0。名前の文字列を示している。
        gr_list = [distin_dict[0], distin_dict[1]]
        log.info("gr_list {}.".format(gr_list))

        # At first, we check difference of genotype between two sample
        # that described at the beginning of each group
        top_smpl_list = [
            glv.conf.g_members_dict[gr_list[0]][0],
            glv.conf.g_members_dict[gr_list[1]][0]]
        log.info("top_smpl_list {}.".format(top_smpl_list))

        # ================================================================
        start = time.time()
        # write out to file
        out_txt_file = distin_dict['variant']['out_path']
        utl.save_to_tmpfile(out_txt_file)

        # ここがparallele化できるか
        # f.writeの最後のflash必要か。
        with open(out_txt_file, mode='a') as f:

            # write header
            f.write("{}\n".format(distin_dict['variant']['hdr_text']))

            # access to vcf using iterater
            for record in vcf_ittr:

                # 1. Skip same GT between top two sample
                if self._skip_same_GT_between_top2sample(
                    record, top_smpl_list) > 0:
                    continue

                # 2. Check GT in your own group
                if self._skip_different_GT_in_own_group(
                    record, top_smpl_list, gr_list) > 0:
                    continue

                # 3. Select different allele combination among 2x2 allele
                asel = AlleleSelect()
                asel.select_diff_allele(record, top_smpl_list, gr_list)

                # skip if pick_mode is different
#                if utl.is_my_pick_mode(
#                    asel.var_type, distin_dict['pick_mode']) != True:
#                    continue
                
                # 4. Save variant information as text file
                for var_type, line in zip(asel.var_types, asel.lines):
                    if utl.is_my_pick_mode(
                        var_type, distin_dict['pick_mode']) == True:
                        f.write("{}\n".format(line))

        log.info("variant {} {}".format(
            utl.elapsed_time(time.time(), start),
            distin_dict['variant']['base_nam']))


    def _skip_different_GT_in_own_group(self, record, tsl, gr_list):

        skip = glv.SKIP_DONT_SKIP

        # check twice, group0, and group1
        for gr_no in range(2):
            # pick sample name belong to a group
            for (sample_no, sample_name) in enumerate(
                glv.conf.g_members_dict[gr_list[gr_no]]):

                if sample_no == 0:
                    continue    # self

                sample0 = tsl[gr_no]
                sample1 = sample_name
                # もし、サンプル間でvariantが見つかった場合は、
                s0_0, s0_1, s1_0, s1_1 = \
                    AlleleSelect.record_call_for_sample(
                        record, sample0, sample1)

                # compare alleles with first sample
                if s0_0 == s1_0 and s0_1 == s1_1:

                    #log.debug("SKIP_SAME_HOMO {},({}){} {}{}/{}{}".format(
                    #    gr_list[gr_no],
                    #    sample_no, sample_name,
                    #    record.call_for_sample[tsl[gr_no]].gt_alleles[0],
                    #    record.call_for_sample[tsl[gr_no]].gt_alleles[1],
                    #    record.call_for_sample[sample_name].gt_alleles[0],
                    #    record.call_for_sample[sample_name].gt_alleles[1]))
                    pass

                else:
                    skip = glv.SKIP_DIFF_INGROUP
                    #log.debug("SKIP_SAME_HOMO {},({}){} {}{}/{}{}".format(
                    #    gr_list[gr_no],
                    #    sample_no, sample_name,
                    #    record.call_for_sample[tsl[gr_no]].gt_alleles[0],
                    #    record.call_for_sample[tsl[gr_no]].gt_alleles[1],
                    #    record.call_for_sample[sample_name].gt_alleles[0],
                    #    record.call_for_sample[sample_name].gt_alleles[1]))
                    return skip

        return skip


    def _skip_same_GT_between_top2sample(self, record, tsl):

        # for REF 20200708
        sample0 = tsl[0]
        sample1 = tsl[1]


        s0_0, s0_1, s1_0, s1_1 = \
            AlleleSelect.record_call_for_sample(record, sample0, sample1)

        skip = glv.SKIP_DONT_SKIP

        # same homo: AA,AA
        if Variant.is_same_homo(s0_0, s0_1, s1_0, s1_1):
            skip = glv.SKIP_SAME_HOMO
            #log.debug("SKIP_SAME_HOMO {}{}/{}{}".format(s0_0,s0_1,s1_0,s1_1))
            return skip

        # same hetero: AB,AB
        if Variant.is_same_hetero(s0_0, s0_1, s1_0, s1_1):
            skip = glv.SKIP_SAME_HETERO
            #log.debug("SKIP_SAME_HETERO {}{}/{}{}".format(
            #    s0_0,s0_1,s1_0,s1_1))
            return skip

        # ./.
        if Variant.is_None(s0_0, s0_1, s1_0, s1_1):
            skip = glv.SKIP_None
            #log.debug("SKIP_None {}{}/{}{}".format(s0_0,s0_1,s1_0,s1_1))
            return skip

        return skip


    @classmethod
    def is_same_gt(cls, s0_0, s0_1, s1_0, s1_1):

        same_gt = False

        # same homo: AA,AA
        if Variant.is_same_homo(s0_0, s0_1, s1_0, s1_1):
            same_gt = True

        # same hetero: AB,AB
        elif Variant.is_same_hetero(s0_0, s0_1, s1_0, s1_1):
            same_gt = True

        return same_gt

    @classmethod
    def is_None(cls, s0_0, s0_1, s1_0, s1_1):
        return s0_0 == None or s0_1 == None or \
               s1_0 == None or s1_1 == None


    @classmethod
    def is_homo_homo(cls, s0_0, s0_1, s1_0, s1_1):
               #  1 == 2           3 == 4
        return s0_0 == s0_1 and s1_0 == s1_1


    @classmethod
    def is_homo_hetero(cls, s0_0, s0_1, s1_0, s1_1):
               #  1 == 2           3 != 4
        return s0_0 == s0_1 and s1_0 != s1_1 or \
               s0_0 != s0_1 and s1_0 == s1_1
               #  1 != 2           3 == 4


    @classmethod
    def is_hetero_hetero(cls, s0_0, s0_1, s1_0, s1_1):
               #  1 != 2           3 != 4
        return s0_0 != s0_1 and s1_0 != s1_1


    @classmethod
    def is_same_homo(cls, s0_0, s0_1, s1_0, s1_1):
        return Variant.is_homo_homo(s0_0, s0_1, s1_0, s1_1) and \
               s0_0 == s1_0
               #  1 == 3


    @classmethod
    def is_same_hetero(cls, s0_0, s0_1, s1_0, s1_1):
        return Variant.is_hetero_hetero(s0_0, s0_1, s1_0, s1_1) and \
               s0_0 == s1_0 and s0_1 == s1_1
               #  1 == 3           2 == 4

    @classmethod
    def is_share(cls, s0_0, s0_1, s1_0, s1_1):
               #  1 == 3          1 == 4
        return s0_0 == s1_0 or s0_0 == s1_1 or \
               s0_1 == s1_0 or s0_1 == s1_1
               #  2 == 3          2 == 4


