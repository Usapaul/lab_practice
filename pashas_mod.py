#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from astropy import units as u
from astropy.modeling.blackbody import blackbody_lambda


# blackbody_lambda(wavelength, temperature)

def Planck_fun(wavelength, temperature):
    return blackbody_lambda(wavelength, temperature).value

