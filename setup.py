#!/usr/bin/env python
# Copyright (c) 2015, Di Zhang <hang.dai@bcm.edu>, Gao Wang <wangow@gmail.com>

# $File: setup.py $
# $LastChangedDate:  $
# $Rev:  $
# copied from gaow's setup file

from distutils.core import setup
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

from src import VERSION
NAME = "VMT"

setup(name = NAME,
    version = VERSION,
    description = "A tool to identify causal variants for Mendelian diseases using sequence data",
    author = "Hang Dai",
    packages = [NAME],
    scripts = ['src/vmt'],
    package_dir = {NAME:'src'},
    cmdclass = {'build_py': build_py },
    data_files=[('format_files', ['format_files/vcf.fmt', 'format_files/VMT_annotation.fmt'])]
)
