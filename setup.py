#/usr/bin/env python

import sys
import os
import shutil

CLASSIFIERS = """
Development Status :: 1 - Planning
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Programming Language :: C
Programming Language :: Python :: 2
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

PKGNAME = "pycoevolve"
DESCR = "A simple python module for coevolutionary metrics."
AUTHOR = "Jacob Mick"
EMAIL = "jam7w2@mail.missouri.edu"
URL = "https://github.com/JakeMick/pycoevolve"
LICENSE = "BSD"
VERSION = "planning-git"


def configuration(parent_package="",top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    config.add_subpackage("pycoevolve")
    config.add_data_files(("pycoevolve/toy_data","*.fa"))
    return config

def setup_package():
    from numpy.distutils.core import setup

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    src_path = local_path
    if sys.version_info[0] == 3:
        src_path = os.path.join(local_path, "build", "py3k")
        sys.path.insert(0, os.path.join(local_path, "tools"))
        import py3tool
        print("Converting to Python3 via 2to3.")
        py3tool.sync_2to3("pycoevolve", os.path.join(src_path, "pycoevolve"))

        site_cfg = os.path.join(local_path, "site.cfg")
        if os.path.isfile(site_cfg):
            shutil.copy(site_cfg, src_path)

    os.chdir(local_path)
    sys.path.insert(0, local_path)
    sys.path.insert(0, os.path.join(local_path, "pycoevolve"))

    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    try:
        setup(name = PKGNAME,
            version = VERSION,
            description = DESCR,
            author = AUTHOR,
            author_email = EMAIL,
            license = LICENSE,
            url = URL,
            classifiers = [row for row in CLASSIFIERS.split("\n") if row],
            configuration = configuration)

    finally:
        del sys.path[0]
        os.chdir(old_path)

if __name__ == "__main__":
    setup_package()
