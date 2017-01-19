# coding=utf-8
#*********************************************************************
# PROGRAM NAME: Intensity.py                                         *
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
# Intensity.py is a class which methods will be used as functions to 
# represent the different intensity distributions at the entrance pla-
# ne.
#
# 'SPATIAL COHERENCE WAVELETS' is a project which focus on unidimensi-
# onal diffractional processes under the nonparaxial modelling of op-
# tical fields. [1]. The simulation is based on the theory of spatial 
# coherence wavelets which purpose is the modelling of optical fields
# in any state of spatial coherence.                                
#*********************************************************************
# METHODS:   
#    slide             : Single slide distribution  
#    double_slide      : Double slide distribution
#    pinhole           : Pin hole distrubution
#    double_pinhole    : Double pinhole distribution
#    sine              : Sine distribution 
#    ronchi_slide      : Ronchi slide distribution
#    straight_line     : DC level distribution with slope variation
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


def slide(x, so, width):
    """
    # NOTE: All lenght units must be in micrometers.
    :param x:            Variable along the domain of the function.
    :param so:           Maximum intensity.
    :param width:        Slide width.
    :return: intensity:  Intensity profile
    """
    ng = 4
    arg = 2*x/width
    support = np.exp(-np.power(arg, ng))
    intensity = so*support
    return intensity



def double_slide(x, so_1, so_2, width_1, width_2, slide_space):
    """
    # NOTE: All lenght units must be in micrometers.
    :param x:            Variable along the domain of the function.
    :param so_1:         Maximum intensity left slide.
    :param so_2:         Maximum intensity right slide.
    :param width_1:      Left slide width.
    :param width_2:      Right slide width.
    :param slide_space:  Space between slides.
    :return: intensity:  Intensity profile
    """
    ng = 1000

    arg_1 = 2*(x + (slide_space + width_1)/2)/width_1
    arg_2 = 2*(x - (slide_space + width_2)/2)/width_2

    slide_1 = np.exp(-np.power(arg_1, ng))
    slide_2 = np.exp(-np.power(arg_2, ng))

    intensity = so_1*slide_1 + so_2*slide_2
    return intensity


def pinhole(x, so, width):
    """
    # NOTE: All lenght units must be in micrometers.
    :param x:            Variable along the domain of the function.
    :param so:           Maximum intensity.
    :param width:        Slide width.
    :return: intensity:  Intensity profile
    """
    ng = 2
    arg = 2*x/width
    support = np.exp(-np.power(arg, ng))
    intensity = so*support
    return intensity


def double_pinhole(x, so_1, so_2, width_1, width_2, slide_space):
    """
    # NOTE: All lenght units must be in micrometers.
    :param x:            Variable along the domain of the function.
    :param so_1:         Maximum intensity.
    :param so_2:         Maximum intensity.
    :param width_1:      Left slide width.
    :param width_2:      Right slide width.
    :param slide_space:  Space between slides.
    :return: intensity:  Intensity profile
    """
    ng = 2

    arg_1 = 2*(x + (slide_space + width_1)/2)/width_1
    arg_2 = 2*(x - (slide_space + width_2)/2)/width_2

    slide_1 = np.exp(-np.power(arg_1, ng))
    slide_2 = np.exp(-np.power(arg_2, ng))

    intensity = so_1*slide_1 + so_2*slide_2
    return intensity


def sine(x, so, width, a, w):
    """
    # NOTE: All lenght units must be in micrometers.
    :param x:            Variable along the domain of the function.
    :param so:           DC level of intensity.
    :param width:        Aperture window width.
    :param a:            Amplitude.
    :param w:            Spatial frequency.
    :return: intensity:  Intensity profile
    """
    wave_length = width/w
    if so >= a and np.mod(w, 2) == 1:
        intensity = so + a*np.sin(2*np.pi*x/wave_length + 0.5*np.pi)
    if so >= a and np.mod(w, 2) == 0:
        intensity = so - a*np.sin(2*np.pi*x/wave_length + 0.5*np.pi)
    else:
        raise ValueError('so >= a, and w integer must be satisfied')
    return intensity


def ronchi_slide(x, so, width, w):
    """
    # NOTE: All lenght units must be in micrometers.
    :param x:            Variable along the domain of the function.
    :param so:           Maximum intensity.
    :param width:        Aperture window width.
    :param w:            Spatial frequency.
    :return: intensity:  Intensity profile
    """
    wave_length = width/w
    if np.mod(w, 2) == 1:
        steps = np.sin(2*np.pi*x/wave_length + 0.5*np.pi) > 0
    if np.mod(w, 2) == 0:
        steps = -np.sin(2*np.pi*x/wave_length + 0.5*np.pi) > 0
    else:
        raise ValueError('Error: w must be an integer')

    intensity = so*steps
    return intensity


def straight_line(x, so_1, so_2, width):
    """
    # NOTE: All lenght units must be in micrometers.
    :param x:            Variable along the domain of the function.
    :param so_1:         Intensity at the left border.
    :param so_2:         Intensity at the right border.
    :param width:        Aperture width.
    :return: intensity:  Intensity profile
    """
    m = (so_2 - so_1)/width
    intensity = m*x + so_1
    return intensity

#*********************************************************************
