# encoding: utf-8
from __future__ import division, print_function, unicode_literals
import os
import numpy as np


def rmse(y, yobs):
    if yobs.size == 0:
        return np.nan
    if np.any(np.isnan(yobs)):
        return np.nan
    return np.sqrt(np.sum((y - yobs)**2) / yobs.size)


def rrmse(y, yobs):
    if yobs.size == 0:
        return np.nan
    if np.any(np.isnan(yobs)):
        return np.nan
    res = y - yobs
    return np.sqrt(np.sum((res/yobs)**2) / yobs.size)


def nrmse(y, yobs):
    if yobs.size == 0:
        return np.nan
    if np.any(np.isnan(yobs)):
        return np.nan
    res = y - yobs
    return np.sqrt(np.sum(res**2) / yobs.size) / np.mean(yobs)

def nse(y, yobs):
    if yobs.size == 0:
        return np.nan
    if np.any(np.isnan(yobs)):
        return np.nan
    
    return 1.0 - np.sum((y - yobs)**2) / np.sum((yobs - np.mean(yobs))**2)
