# coding=utf-8
#*********************************************************************
# PROGRAM NAME: ToRun.py                                             *
#                                                                    *
#*********************************************************************
# AUTHOR: David Alejandro Vargas O                                   *
#         david.vargas@geophysik.uni-muenchen.de                     *
# AFFILIATION: LMU Ludwig Maximilian University of Munich            *
# DATE: 27.03.2015                                                   *
#*********************************************************************
# DESCRIPTION:
# This file is a part of the  Project SPATIAL COHERENCE WAVELETS.

# ToRun is the principal file of the project SPATIAL COHERENCE WAVELE-
# TS. here, the methods that compute the marginal power spectrum, the
# Intensity distribution function, and the degree of spatial coheren-
# ce function are called. these functions computed under specific phy-
# sical and geometric parameters became the key point in the analysis
# of optical fields in arbitrary states of spacial coherence that pro-
# pagate between two planes, aperture or entrance plane 'AP' and obser
# vation or exit plane 'OP' according with the theory of spatial cohe-
# rence wavelets [1].

# 'SPATIAL COHERENCE WAVELETS' is a project which focus on unidimensi-
# onal diffractional processes under the nonparaxial modelling of op-
# tical fields. [1]. The simulation is based on the theory of spatial 
# coherence wavelets which purpose is the modelling of optical fields
# in any state of spatial coherence.                               
#*********************************************************************
# PARAMETERS:                                                        
#    so          :  Maximum intensity.
#    wide        :  Slide wide.
#    so_1        :  Maximum intensity left slide.
#    so_2        :  Maximum intensity right slide.
#    wide_1      :  Left slide wide.
#    wide_2      :  Right slide wide.
#    slide_space : Space between slides.
#    a           : Amplitude.
#    w           : Spatial frequency.
#
#    sources       : Point sources amount.
#    pixel_number  : sampling pixel amount in xa
#    wavelength    : Optical field wavelength.
#    z             : Aperture plane - Exit plane distance.
#    aperture_size : Aperture size
#    exit_size     : Exit size
#
#    sigma   : Gaussian standard deviation.
#    gamma   : Lorentzian standard deviation.
#
# INPUT:   input1, input2                                           
# OUTPUT:  output1, output2
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

# IMPORT MODULES
import Intensity
import CoheDegree
import MarPowSpec

import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------------------------------------------
# SET PARAMETERS

# COORDINATES AXIS: xi_A(AP APERTURE PLANE), r_A(OP OBSERVATION PLANE)
wavelength = 0.632     # Optical field wavelength um.
z = 50                 # AP-OP distance um. 
sources = 20           # Point sources amount in xi_A.
pixel_number = 1024    # sampling pixel amount in r_A.

alpha = np.sqrt(0.1)  # consider 99% of the free space diffraction envelope
AP_width = 10          # AP Aperture plane width um.
OP_width = np.sqrt(np.power((1+np.sqrt(1+8*alpha))/(4*alpha),2)-1)*z + AP_width 
xi_A = np.linspace(-0.5 * AP_width, 0.5 * AP_width, 2*sources - 1)  # AP Axis
r_A = np.linspace(-0.5 * OP_width, 0.5 * OP_width, pixel_number)    # OP Axis

# INTENSITY
so = 1             # Maximum intensity.
#width = 10        # Slide width.
#so_1 = 100        # Maximum intensity left slide.
#so_2 = 100        # Maximum intensity right slide.
#width_1 = 40      # Left slide width.
#width_2 = 40      # Right slide width.
#slide_space = 60  # Space between slides.
#a = 100           # Amplitude.
#w = 4             # Spatial frequency.

# SET COHERENCE DEGREE PARAMETERS
sigma = 1000000    # Gaussian standard deviation.
#gamma = 1000      # Lorentzian standard deviation.

##FBeta = 5
#wavelength = wide/((sources-1)*FBeta);

#---------------------------------------------------------------------------
# SET INTENSITY
intensity = lambda x: Intensity.slide(x, so, AP_width)

# SET COHERENCE
cohe_degree = lambda x: CoheDegree.gauss_cohe(x, sigma)

# MARGINAL POWER SPECTRUM
MPSpectrum = lambda xi_A, r_A, z: MarPowSpec.nonparaxial(xi_A, r_A, z,
                                   wavelength, intensity, cohe_degree)
        
#---------------------------------------------------------------------------
# ANALYSIS OF LIGHT PROPAGATION                                  
# COMPUTE MARGINAL POWER SPECTRUM
MPS_rad, MPS_vir = MPSpectrum(xi_A, r_A, z)
Mar_Pow_Spectrum = MPS_rad + MPS_vir

# COMPUTE POWER SPECTRUMS 
PS_rad_OP = MPS_rad.sum(axis=1)  # Radiant power spectrum at observation plane 
PS_vir_OP = MPS_vir.sum(axis=1)  # Virtual power spectrum at observation plane 
PS_OP = PS_rad_OP + PS_vir_OP    # Power spectrum at observation plane
PS_AP = Mar_Pow_Spectrum.sum(axis=0) # Power spectrum at aperture plane AP

# COMPUTE NORMALIZED POWER SPECTRUMS
N_PS_rad_OP = PS_rad_OP/np.max(PS_rad_OP) # PS_rad_OP Normalized
N_PS_vir_OP = PS_vir_OP/np.max(PS_vir_OP) # PS_vir_OP Normalized
N_PS_OP = PS_OP/np.max(PS_OP) # PS_OP Normalized
N_PS_AP = PS_AP/np.max(PS_AP) # PS_AP Normalized

#---------------------------------------------------------------------------        
# SET GRAPHICS
#MARGINAL POWER SPECTRUM, ENTRANCE AND EXIT POWER SPECTRUM
plt.figure(figsize=(16,10), dpi=80) 
plt.subplot(1, 2, 1)
plt.imshow(Mar_Pow_Spectrum, cmap='gray', aspect='auto',
           extent=[-0.5*AP_width, 0.5*AP_width,
                   -0.5*OP_width, 0.5*OP_width])
plt.title('Marginal Power Spectrum')
plt.xlabel(r'$ \xi_A \ (\mu m)$')
plt.ylabel(r'$ r_A \ (\mu m)$')

plt.subplot(1, 2, 2)
#plt.plot(xa, s_xaNormal)
plt.fill_between(r_A, 0, N_PS_OP) # 'r', alpha=0.7
plt.plot(r_A, N_PS_rad_OP, 'r--')
plt.title('Power Spectrum - Observation plane')
plt.xlabel(r'$ r_A \ (\mu m) $')
plt.ylabel(r'$ S(r_A) $')

#---------------------------------------------------------------------------
## ANALYSIS OF LIGHT PROPAGATION                                  
#longitude = np.linspace(0.1, 1000, 1000)
#power = []
#for z in longitude:
#    # COMPUTE MARGINAL POWER SPECTRUM
#    mpeReal, mpeVirtual = MarPowSpec(xia, xa, z)
#    #Mar_Pow_Spectrum = mpeReal + mpeVirtual
#
#    # POWER SPECTRUM
#    s_xaReal = mpeReal.sum(axis=1) # REAL LAYER POWER SPECTRUM
#    s_xaRealNormal = s_xaReal/np.amax(s_xaReal);
#    #s_xaVirtual = mpeVirtual.sum(axis=1) # VIRTUAL LAYER POWER SPECTRUM
#    #s_xa = s_xaReal + s_xaVirtual # EXIT POWER SPECTRUM 
#    #s_xia = Mar_Pow_Spectrum.sum(axis=0) # ENTRANCE POWER SPECTRUM
#    #s_xaNormal=s_xa/np.amax(s_xa);
#    power.append(s_xaRealNormal)    
#power = np.asanyarray(power)
#power = power>0.01
##power = power[:, 0:1024] 
#  
#plt.imshow(power.T, cmap='hot', aspect='auto',
#           extent=[0, z, 0, 0.5*exit_size])
#    #title(sprintf('Potencia espectral a lo largo del eje z'))
#    #xlabel('Z  [um]')
#    #ylabel('\xi_A  [um]')                     
#---------------------------------------------------------------------------
## SET GRAPHICS
##MARGINAL POWER SPECTRUM, ENTRANCE AND EXIT POWER SPECTRUM
#plt.figure(figsize=(16,10), dpi=80) 
#plt.subplot(1, 3, 1)
#plt.imshow(Mar_Pow_Spectrum, cmap='gray', aspect='auto',
#           extent=[-0.5*aperture_size, 0.5*aperture_size,
#                   -0.5*exit_size, 0.5*exit_size])
#plt.title('Marginal Power Spectrum')
#plt.xlabel(r'$ \xi_A \ (\mu m)$')
#plt.ylabel(r'$ r_A \ (\mu m)$')
#
#plt.subplot(1, 3, 2)
#plt.plot(xa, s_xaNormal)
##plt.fill_between(xa, 0, s_xaNormal) # 'r', alpha=0.7
#plt.plot(xa, s_xaRealNormal)
#plt.title('Exit Power Spectrum')
#plt.xlabel(r'$ r_A \ (\mu m) $')
#plt.ylabel(r'$ S(r_A) $')
#
#plt.subplot(1, 3, 3)
#plt.plot(xa, s_xaReal, xa, s_xaVirtual)
#plt.plot(xa, -s_xaReal, 'b')
#plt.title('Entrance Power Spectrum')
#plt.xlabel(r'$ \xi_A \ (\mu m)$')
#plt.ylabel(r'$ S(\xi_A) $')
#plt.show()
#*********************************************************************

