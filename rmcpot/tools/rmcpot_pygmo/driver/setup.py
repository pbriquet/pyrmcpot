#!/usr/bin/env python3

from distutils.core import setup

setup(name='rmcpot_driver',version='0.1',
      requires=['pygmo','mumpy','rmcpot_prob'],
      author='Roberto Veiga',packages=['.'])
