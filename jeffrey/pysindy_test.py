# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 15:51:01 2020

@author: STARS
"""

import numpy as np
import pysindy as ps

t = np.linspace(0, 1, 100)
x = 3 * np.exp(-2 * t)
y = 0.5 * np.exp(t)
X = np.stack((x, y), axis=-1)  # First column is x, second is y

model = ps.SINDy(feature_names=["x", "y"])
model.fit(X, t=t)

model.print()