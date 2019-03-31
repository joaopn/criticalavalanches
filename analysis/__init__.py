# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-24 19:03:08
# @Last Modified by:   joaopn
# @Last Modified time: 2019-03-31 18:55:56

from analysis import avalanche, plot

import numpy as np
import os
import matplotlib
if os.environ.get('DISPLAY', '') == '':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
