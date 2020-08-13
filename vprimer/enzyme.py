# -*- coding: utf-8 -*-


# http://rebase.neb.com/rebase/rebase.enz.html
# https://international.neb.com/tools-and-resources/selection-charts/isoschizomers

import sys
import os
import errno

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
log = LogConf.open_log(__name__)

import Bio.Restriction.Restriction_Dictionary as ResDict

class Enzyme(object):

    def __init__(self):

        self.enzyme_list = list()

    def read_enzyme_file(self):

        enz_list = list()

        for enzyme_file in glv.conf.enzyme_files_list:
            log.info("{}".format(enzyme_file))
            with open(enzyme_file, "r") as f:
                # iterator
                for r_liner in f:
                    r_line = r_liner.strip()    # cr, ws

                    if r_line.startswith('#') or r_line == '':
                        continue
                    r_line = utl.strip_hash_comment(r_line)

                    if r_line not in ResDict.rest_dict:
                        log.critical("your ENZYME {} is not in list.".format(
                            r_line))
                        sys.exit(1)

                    enz_list.append(r_line)

        self.enzyme_list = list(set(enz_list))
        self.enzyme_list.sort(key=None, reverse=False)
#        log.info("num={} {}".format(
#            len(self.enzyme_list), self.enzyme_list))

    @classmethod
    def prepare_enzyme(cls):

        #enzyme_files=
        #    data/DaiichiSeiyaku_NEB_07.txt;\
        #    data/RE_labo_v1.0_20150403.txt;
        #    data/Takara07.txt

        enzyme_files_list = \
            glv.conf.enzyme_files_str.split(';')

        log.info("enzyme_files_list {}".format(enzyme_files_list))

        enzyme_file = ''

        for enzyme_file_user in enzyme_files_list:

            # originally absolute path
            if enzyme_file_user.startswith('/'):
                enzyme_file = enzyme_file_user
            else:
                enzyme_file = "{}/{}".format(glv.conf.cwd, enzyme_file_user)

            log.info("enzyme_file {}".format(enzyme_file))

            # slink
            basename_user = os.path.basename(enzyme_file_user)
            enzyme_file_slink_system = "{}/{}".format(
                glv.conf.ref_dir, basename_user)
            # save as conf global
            glv.conf.enzyme_files_list.append(enzyme_file_slink_system)

            if os.path.isfile(enzyme_file_slink_system):
                log.info("{} exist.".format(enzyme_file_slink_system))
            else:
                log.info("os.symlink {} {}.".format(
                    enzyme_file, enzyme_file_slink_system))

                os.symlink(enzyme_file, enzyme_file_slink_system)

        log.info("enzyme_files {}".format(glv.conf.enzyme_files_list))


