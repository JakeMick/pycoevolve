#/usr/bin/env python
import os
from distutils.core import setup

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

setup(name=PKGNAME,
      version=VERSION,
      description=DESCR,
      author=AUTHOR,
      author_email=EMAIL,
      url=URL,
      packages=["pycoevolve", "pycoevolve.io",
          "pycoevolve.toy_data","pycoevolve.utils"]
     )

