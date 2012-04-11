#/usr/bin/env python
from distutils.core import setup

CLASSIFIERS = """
Development Status :: 1 - Planning
Intended Audience :: Science/Research
License :: OSI Approved :: Simplified BSD
Programming Language :: C
Programming Language :: Python :: 2
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
LICENSE = "Simplified BSD"
VERSION = "planning-git"

setup(  name = PKGNAME,
        version = VERSION,
        description = DESCR,
        author = AUTHOR,
        author_email = EMAIL,
        url = URL,
        classifiers = [cl for cl in CLASSIFIERS.split("\n")
            if cl is not ""],
        packages=["pycoevolve",
            "pycoevolve.seqs",
            "pycoevolve.seqs.tests",
            "pycoevolve.toy_data",
            "pycoevolve.toy_data.tests",
            "pycoevolve.utils",
            "pycoevolve.metrics"],
        package_data = {"pycoevolve" :
            ["toy_data/toy_fastas/*fa"]}
        )

