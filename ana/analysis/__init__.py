# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-24 19:03:08
# @Last Modified by:   joaopn
# @Last Modified time: 2019-04-01 02:05:54

from analysis import avalanche, plot, fitting

import numpy as np
import os
import matplotlib
if os.environ.get('DISPLAY', '') == '':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
