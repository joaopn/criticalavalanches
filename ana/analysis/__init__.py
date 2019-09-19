# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-24 19:03:08
# @Last Modified by:   Joao
# @Last Modified time: 2019-07-06 10:47:53


import numpy as np
import os
import matplotlib
if os.environ.get('DISPLAY', '') == '':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

from analysis import avalanche, plot, fitting, parser
from analysis.dataset import *

