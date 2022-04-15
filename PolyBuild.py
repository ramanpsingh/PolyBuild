#!/usr/bin/python

import sys
import math
from math import *
import argparse

# Argument Parser
parser = argparse.ArgumentParser(description='Create a custom Pluronic coordinate and topology file')
parser.add_argument( "-peo", "--peo", type=int, default=99, help='Number of PEO beads in one block')
parser.add_argument( "-ppo", "--ppo", type=int, default=69, help='Number of PPO beads in the central block')
parser.add_argument( "-bl", "--bondlength", type=float, default=0.280, help='Bond length')
parser.add_argument( "-o", "--output",   type=str,   default=None,   help='Output file name, default: generate automatically' )
args = parser.parse_args()

eo = args.peo
po = args.ppo
l = args.bondlength

block1=eo
block2=eo+po
block3=eo+po+eo

if args.output == None:
	output = "Pluronic-EO"+str(eo)+"-PO"+str(po)+"-EO"+str(eo)
else:
        output = args.output

numatoms = block3

# ang1 and ang2 are angles (in degress) between EO-EO and PO-PO beads.
ang1=155
ang2=144

print( "Generated Martini model for Pluronic with EO%s-PO%s-EO%s" % (eo, po, eo))

# CREATE GRO FILE

# Opens the gro file for writing

structure_file = open(output+".gro", 'w')


# Writes header for gro file

structure_file.write( "Pluronic-EO%s-PO%s-EO%s\n" % (eo, po, eo) )
structure_file.write( "  %3d\n" % (numatoms) )


# Writes coordinates for beads
# The polymer chain is aligned along the x-axis.

z = 0

for i in range(1,block1+1):
	x_eo=i*l*sin(ang1/2)
	if i%2==1:
		y_eo=0
	else:
		y_eo=l*cos(ang1/2)
	structure_file.write(  "%5d%-5s E%03d%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "  PLU", i, i, x_eo, y_eo, z, 0, 0, 0) )	

for i in range(block1+1,block2+1):
	x_po=(i-eo)*l*sin(ang2/2) + (eo*l*sin(ang1/2))
	if i%2==1:
		y_po=0
	else:
		y_po=l*cos(ang2/2)
	structure_file.write(  "%5d%-5s P%03d%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "  PLU", i, i, x_po, y_po, z, 0, 0, 0) )	

for i in range(block2+1,block3+1):
	x_eo=(i-eo-po)*l*sin(ang1/2) + (eo*l*sin(ang1/2)) + (po*l*sin(ang2/2))
	if i%2==1:
		y_eo=0
	else:
		y_eo=l*cos(ang1/2)
	structure_file.write(  "%5d%-5s E%03d%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "  PLU", i, i, x_eo, y_eo, z, 0, 0, 0) )	

# Box dimensions

xdim = 2*eo*l*sin(ang1/2) + po*l*sin(ang2/2) + 2
ydim = round(l*cos(ang1))+2
zdim = 2

structure_file.write(  "   %u  %u  %u\n" % (xdim, ydim, zdim) )

structure_file.close()


# CREATE TOPOLOGY FILE

# Opens the itp file for writing

topology_file = open(output+".itp", 'w')

# Writes header for the itp file

topology_file.write( "; \n; Martini topology for Pluronic with EO%s-PO%s-EO%s\n" % (eo, po, eo) )
topology_file.write( "; \n; Topology generated using PolyBuild.py \n" )
topology_file.write( "; Written by Raman Preet Singh \n")
topology_file.write( "; \n; \n; The bead types and bonded parameters are based on: \n")
topology_file.write( "; G. Perez-Sanchez, F.A. Vicente, N. Schaeffer, I.S. Cardoso, S.P.M. Ventura, M. Jorge, J.A.P. Coutinho \n")
topology_file.write( "; Rationalizing the Phase Behavior of Triblock Copolymers through Experiments and Molecular Simulations \n")
topology_file.write( "; J. Phys. Chem. C 2019, 123, 21224−21236\n")
topology_file.write( "; DOI: 10.1021/acs.jpcc.9b04099 \n;\n\n")
topology_file.write( "[ moleculetype ]\n" )
topology_file.write( "; Name	 nrexcl\n" )
topology_file.write( "PLU  1\n")


# Atoms

topology_file.write( "\n[ atoms ]\n" )
topology_file.write( "; id	 type	 resnr	 residue	 atom	 cgnr	 charge	 mass\n" )
for i in range(1,block1+1):
	beadtype="SP1"
	topology_file.write( "%3d    %4s   1   PLU    E%03d     %3d\n" % (i, beadtype, i, i) )
for i in range(block1+1,block2+1):
	beadtype="SC3"
	topology_file.write( "%3d    %4s   1   PLU    P%03d     %3d\n" % (i, beadtype, i, i) )
for i in range(block2+1,block3+1):
	beadtype="SP1"
	topology_file.write( "%3d    %4s   1   PLU    E%03d     %3d\n" % (i, beadtype, i, i) )

# Bonds

topology_file.write( "\n[ bonds ]\n" )
topology_file.write( "; i	 j	  funct	 length	 force\n" )

for i in range(1,block1):
	m = i
	n = i + 1
	topology_file.write( "     %3d     %3d       1   %4.3f    8000\n" % (m, n, l) )
for i in range(block1,block2+1):
	m = i
	n = i + 1
	topology_file.write( "     %3d     %3d       1   %4.3f    5000\n" % (m, n, l) )
for i in range(block2+1,block3):
	m = i
	n = i + 1
	topology_file.write( "     %3d     %3d       1   %4.3f    8000\n" % (m, n, l) )

# Angles

topology_file.write( "\n[ angles ]\n" )
topology_file.write( "; i	 j	 k	 funct	 angle	 force\n" )

for i in range(1,block1):
	m = i
	n = i + 1
	o = i + 2
	topology_file.write( "     %3d     %3d     %3d       2     %3.1f     40.0\n" % (m, n, o, ang1) )
for i in range(block1,block2):
	m = i
	n = i + 1
	o = i + 2
	topology_file.write( "     %3d     %3d     %3d       2     %3.1f     40.0\n" % (m, n, o, ang2) )
for i in range(block2,block3-1):
	m = i
	n = i + 1
	o = i + 2
	topology_file.write( "     %3d     %3d     %3d       2     %3.1f     40.0\n" % (m, n, o, ang1) )

