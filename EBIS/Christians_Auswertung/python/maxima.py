#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys

data = np.loadtxt(sys.argv)

for i in range(1, 10):
	print(data[i,1])
