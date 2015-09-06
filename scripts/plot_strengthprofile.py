#!/usr/bin/env python
# plot_strengthprofile.py ---
#
# Plot strength vs. depth for a number of different background temperature
# gradients and activation energies
#
# Author: Iain W. Bailey
# Created: Wed Jun 26 20:44:15 2013 (-0400)
# Last-Updated:
#           By:
#     Update #: 61
#

import sys # libraries used
import numpy as NP
from math import sqrt, log10
import matplotlib.pylab as PLT

# Define global variables
w = 0.1 # Fault width in m
Rg = 8.3144 # gas constant
A = 6.31e-20 # arrhenius amplitude
fs = 0.75 # Friction coefficient
tau0 = 6.0 # Cohesion in MPa

maxz = 17.5; # max depth in km
nz = 101; # number of depths to sample
dsdz = 18.0 # gradient of effective normal stress
Tsurf = 293.0 # surface temperature
vc = 35.0e-3 /(60*60*24*365.25)# velocity in m/s
n = 3.0 # stress exponent
ofile = 'strengthprofile.png'

#-------------------------------------------------------------------------------
def staticStrength( dsdz, fs, z ):
    """
    Calculate the static strength given the cohesion, effective normal
    stress gradient, coefficient of friction and depth
    """
    return tau0 + dsdz*fs*z

#-------------------------------------------------------------------------------
def creepStrength( vc, E, T ):
    """
    Calculate the stress for creep at velocity vc and activation
    energy E
    """
    return 1.0E-6*( (vc/(w*A)) * NP.exp( E/(Rg*T) ))**(1.0/n)

#-------------------------------------------------------------------------------

depths = NP.linspace( 0, maxz, nz)
taus = staticStrength( dsdz, fs, depths )

fig = PLT.figure()
DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*.75, DefaultSize[1]*1.5 ) )
ax = fig.add_subplot(111)
ax.plot( taus, depths )

E = 110.0e3 # Activation energy in J
dTdz = 15.0 # background temperature gradient
temperature = Tsurf + dTdz*depths
tauc = creepStrength( vc, E, temperature )
p2 = ax.plot( tauc, depths, label="E=110|dTdz=15" )

E = 130.0e3 # Activation energy in J
dTdz = 20.0 # background temperature gradient
temperature = Tsurf + dTdz*depths
tauc = creepStrength( vc, E, temperature )
p1 = ax.plot( tauc, depths, label="E=130|dTdz=20" )

E = 130.0e3 # Activation energy in J
dTdz = 25 # background temperature gradient
temperature = Tsurf + dTdz*depths
tauc = creepStrength( vc, E, temperature )
p1 = ax.plot( tauc, depths, label="E=130|dTdz=25" )

E = 135.0e3 # Activation energy in J
dTdz = 20.0 # background temperature gradient
temperature = Tsurf + dTdz*depths
tauc = creepStrength( vc, E, temperature )
p1 = ax.plot( tauc, depths, label="E=135|dTdz=20" )

E = 150.0e3 # Activation energy in J
dTdz = 30.0 # background temperature gradient
temperature = Tsurf + dTdz*depths
tauc = creepStrength( vc, E, temperature )
p3 = ax.plot( tauc, depths, label="E=150|dTdz=30" )

ax.grid(True)
ax.set_ylabel('Depth [km]')
ax.set_xlabel('Strength [MPa]')
ax.axis([0, max(taus), 0.0, maxz] )
ax.invert_yaxis()
ax.legend()
ax.set_aspect(25)

#PLT.show()

PLT.savefig( ofile ,bbox_inches='tight')
sys.stderr.write('Output written to %s\n' % ofile )



# # plot temp profile
# ax = fig1.add_subplot(221)

# #
# ax.plot(Temp1, depths)

# # plot strength profile
# ax = fig1.add_subplot(222)
# ax.plot(strength1, depths)
# ax.invert_yaxis()
# ax.grid(True)





# plot_strengthprofile.py ends here
