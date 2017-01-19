# coding=utf-8
#*********************************************************************
# PROGRAM NAME: MarPowSpec.py                                        *
#                                                                    *
#*********************************************************************
# AUTHOR: David Alejandro Vargas O                                   *
#         david.vargas@geophysik.uni-muenchen.de                     *
# AFFILIATION: LMU Ludwig Maximilian University of Munich            *
# DATE: 27.03.2015                                                   *
#*********************************************************************
# DESCRIPTION:
# This file is a part of the  Project SPATIAL COHERENCE WAVELETS.
# 
# MarPowSpec.py is a class which methods will be used to compute the
# marginal power spectrum corresponding to the propagation of optical
# fields in any state of spatial coherence. There are basically two 
# options, Marginal power spectrums that model progation of optical 
# fields under the aproach of the paraxial aproximation, and those 
# which use the nonparaxial one. This class offer solutions for both 
# cases.
#
# 'SPATIAL COHERENCE WAVELETS' is a project which focus on unidimensi-
# onal diffractional processes under the nonparaxial modelling of op-
# tical fields. [1]. The simulation is based on the theory of spatial 
# coherence wavelets which purpose is the modelling of optical fields
# in any state of spatial coherence.                                
#*********************************************************************
# METHODS:   
#    paraxial        : Marginal power spectrum paraxial approach  
#    nonparaxial     : Margianal power spectrum nonparaxial approach
#    d3paraxial      : 3D Margial power spectrum paraxial approach
#    d3nonparaxial   : 3D Marginal power spectrum nonparaxial approach
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


def paraxial(xi_A, r_A, z, wavelength, intensity, coherence_degree):
    """
    # NOTE: All lenght units must be in micrometers.

    :param xi_A:               Entrance coordinate xi_A.
    :param r_A:                Observation coordinate r_A.
    :param wavelength:         Optical field wavelength.
    :param z:                  Aperture plane - Observation plane distance.
    :param intensity:          Intensity distribution
    :param coherence_degree:   Coherence degree
    :return: real_layer:       Marginal Power Spectrum's real layer
    :return: virtual_layer:    Marginal Power Spectrum's virtual layer
    """

    # BASE MATRIX
    nf = int((np.size(xi_A) + 1)/2)   # Point sources amount.
    n = np.size(r_A)  # number of pixels at the Observation plane OP
    m = 2 * nf - 1  # number of pixels at the aperture plane AP
    grid, nonused = np.meshgrid(np.arange(1, m + 1), np.arange(1, n + 1))

    # COORDINATE MATRICES OF THE APERTURE xi_A AND Observation r_A PLANE
    xi_A_matrix, r_A_matrix = np.meshgrid(xi_A, r_A)

    # WAVELENGTH, DISCRETIZATION OF SLIDE
    k = (2 * np.pi) / wavelength
    AP_width = 2 * xi_A.max()   # AP width 
    pixel_size = AP_width / (2 * nf - 2)

    # COMPUTE REAL LAYER
    real_layer = intensity(xi_A_matrix) * np.mod(grid, 2)

    # COMPUTE VIRTUAL LAYER
    virtual_layer = np.zeros((n, m))
    structure = np.mod(grid, 2) <= 0
    structure = 1*structure    # convert boolean array to int array

    if nf >= 2:
        for family in range(2, nf+1):
            structure[:, 0:family - 1] = 0
            structure[:, (2 * nf - 1) - (family - 2):] = 0
            
            xi_d = 2 * (family - 1) * pixel_size
                             
            rad_1 = np.sqrt(intensity(xi_A_matrix + 0.5*xi_d))
            rad_2 = np.sqrt(intensity(xi_A_matrix - 0.5*xi_d))
            phase = (k/z)* xi_d * xi_A_matrix + (k/z) * xi_d * r_A_matrix
            cosine = np.cos(phase)

            source = 2* coherence_degree(xi_d) * rad_1 * rad_2 * cosine
            
            virtual = structure * source
            virtual_layer = virtual_layer + virtual
    
            structure = structure == 0  # invert zeros for ones.
            # convert boolean array to int array
            structure = 1*structure
    else:
        raise ValueError('np must be mayor or equal 2')
    return real_layer, virtual_layer



def nonparaxial(xi_A, r_A, z, wavelength, intensity, coherence_degree):
    """
    # NOTE: All lenght units must be in micrometers.
    :param xi_A:               Entrance coordinate xi_A.
    :param r_A:                Observation coordinate r_A.
    :param wavelength:         Optical field wavelength.
    :param z:                  Aperture plane - Observation plane distance.
    :param intensity:          Intensity distribution
    :param coherence_degree:   Coherence degree
    :return: real_layer:       Marginal Power Spectrum's real layer
    :return: virtual_layer:    Marginal Power Spectrum's virtual layer
    """

    # BASE MATRIX
    nf = int((np.size(xi_A) + 1)/2)   # Point sources amount.
    n = np.size(r_A)  # number of pixels at the Observation plane OP
    m = 2 * nf - 1  # number of pixels at the aperture plane AP
    grid, nonused = np.meshgrid(np.arange(1, m + 1), np.arange(1, n + 1))
    
    # COORDINATE MATRICES OF THE APERTURE xi_A AND Observation r_A PLANE
    xi_A_matrix, r_A_matrix = np.meshgrid(xi_A, r_A)

    # WAVELENGTH, DISCRETIZATION OF SLIDE
    k = (2 * np.pi) / wavelength
    AP_width = 2 * xi_A.max()   # AP width
    pixel_size = AP_width / (2 * nf - 2)

    # COMPUTE REAL LAYER
    radicand = np.power(z, 2) + np.power((r_A_matrix - xi_A_matrix), 2)
    lorentz_real = np.power((z + np.sqrt(radicand)) / radicand, 2)
    factor1 = intensity(xi_A_matrix) / (4 * np.power(wavelength, 2)) 
    real_layer = factor1 * lorentz_real * np.mod(grid, 2)

    # COMPUTE VIRTUAL LAYER
    virtual_layer = np.zeros((n, m))
    structure = np.mod(grid, 2) <= 0
    structure = 1*structure    # convert boolean array to int array

    if nf >= 2:
        for family in range(2, nf+1):
            structure[:, 0:family - 1] = 0
            structure[:, (2 * nf - 1) - (family - 2):] = 0

            xi_d = 2 * (family - 1) * pixel_size

            term_1 = np.power(z, 2)
            term_2 = np.power(r_A_matrix - xi_A_matrix, 2)
            term_3 = (np.power(xi_d, 2) / 4)
            term = term_1 + term_2 + term_3

            radicand_1 = term - (r_A_matrix - xi_A_matrix) * xi_d
            radicand_2 = term + (r_A_matrix - xi_A_matrix) * xi_d

            lor_1 = (z + np.sqrt(radicand_1)) / radicand_1
            lor_2 = (z + np.sqrt(radicand_2)) / radicand_2

            rad_1 = np.sqrt(intensity(xi_A_matrix + 0.5*xi_d))
            rad_2 = np.sqrt(intensity(xi_A_matrix - 0.5*xi_d))
            
            factor2 = (1 / (2 * np.power(wavelength, 2)))
            lorentz_virtual = factor2 * rad_1 * rad_2 * lor_1 * lor_2
            phase = k * (np.sqrt(radicand_1) - np.sqrt(radicand_2))
            cosine = np.cos(phase)

            source = coherence_degree(xi_d) * lorentz_virtual * cosine

            virtual = structure * source
            virtual_layer = virtual_layer + virtual  

            structure = structure == 0  # invert zeros for ones.
            # convert boolean array to int array
            structure = 1*structure
            
    else:
        raise ValueError('np must be mayor or equal 2')        
    return real_layer, virtual_layer
    
#*********************************************************************
