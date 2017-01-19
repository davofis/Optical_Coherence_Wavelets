# coding=utf-8
#*********************************************************************
# PROGRAM NAME: CoheDegree.py                                        *
#                                                                    *
#*********************************************************************
# AUTHOR: David Alejandro Vargas O                                   *
#         david.vargas@geophysik.uni-muenchen.de                     *
# AFFILIATION: LMU Ludwig Maximilian University of Munich            *
# DATE: 01.04.2015                                                   *
#*********************************************************************
# DESCRIPTION:
# This file is a part of the  Project SPATIAL COHERENCE WAVELETS.
#
# CoheDegree.py is a class which methods will be used as functions
# to represent different degrees of spatial coherence distributions
# at the entrance plane.
#
# 'SPATIAL COHERENCE WAVELETS' is a project which focus on unidimensi-
# onal diffractional processes under the nonparaxial modelling of op-
# tical fields. [1]. The simulation is based on the theory of spatial 
# coherence wavelets which purpose is the modelling of optical fields
# in any state of spatial coherence.                                
#*********************************************************************
# METHODS:   
#    gauss_cohe     : Gaussian spacial coherence degree 
#    lorentz_cohe   : Lorentzian spacial coherence degree 
#
# INPUT: Parameters requiered for each method are expecified indivi-
#        dually inside each method.    
# OUTPUT: Same as above.
#
# NOTE: All lenght units must be in micrometers.
#                                         
# REFERENCES: ********************************************************
# [1] Castañeda R, Sucerquia J. Non-approximated numerical modeling  *
# of propagation of light in any state of spatial coherence          *
# [2] R. Castañeda,* D. Vargas, E. Franco. Discreteness of the set of*
# radiant point sources: a physical feature of the second-order wave *
# -fronts                                                            *
# [3] R. Castañeda,* D. Vargas, E. Franco. Spatial coherence of light*
# and a fundamental discontinuity of classical second-order wave     *
# fronts                                                             *
#*********************************************************************
#     COPYRIGHT                                                      *
# Copyright (C) 2015  David Alejandro Vargas Otalora                 *
#                                                                    *
# This program is free software: you can redistribute it and/or modif*
# y it under the terms of the GNU General Public License as published*
# by the Free Software Foundation, either version 3 of the License,  *
# or(at your option) any later version.                              *
#                                                                    *
# This program is distributed in the hope that it will be useful,    *
# but WITHOUT ANY WARRANTY; without even the implied warranty of     *
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      *
# GNU General Public License for more details.                       *
#                                                                    *
# You should have received a copy of the GNU General Public License  *
# along with this program.  If not, see http://www.gnu.org/licenses/ *
#*********************************************************************

# IMPORT LIBRARIES
import numpy as np


def gauss_coherence(x, sigma):
    """
    # NOTE: All lenght units must be in micrometers.
    :param x:            Variable along the domain of the function.
    :param sigma:        Standard deviation.
    :return: gauss_coh:  coherence profile
    """
    d_pow = np.power(x, 2)
    sigma_pow = 2 * np.power(sigma, 2)
    gauss_coherence = np.exp(- d_pow / sigma_pow)
    return gauss_coherence


def lorentz_coherence(x, gamma):
    """
    # NOTE: All lenght units must be in micrometers.
    :param x:            Variable along the domain of the function.
    :param gamma:        Standard deviation.
    :return: lorentz_coh:  coherence profile
    """
    d_pow = np.power(x, 2)
    gamma_pow = np.power(gamma, 2)
    lorentz_coherence = gamma_pow / (gamma_pow + d_pow)
    return lorentz_coherence

#*********************************************************************
