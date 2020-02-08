# Homework 3 Solutions
# Computing the Mass of the Local Group
# G. Besla


# import relevant modules
import numpy as np
import astropy.units as u
from ReadFile import Read



def ComponentMass(filename,PType):
# function to compute the mass of all particles of a given type for a given galaxy 
# input:  filename (str) ,  Particle type (1,2,3)
# output: Mass in units of 1e12 Msun. 
    
     # read in the file 
    time, total, data = Read(filename)
    
    #create an array to store indexes of particles of desired Ptype
    index = np.where(data['type'] == PType)
    
    # create new arrays with the m of just the Ptype desired.
    mnew = data['m'][index]*1e10*u.Msun
  
    
    # sum the mass of all particles to get total
    # round to 3 decimal places
    mass = np.round(sum(mnew)/1e12,3)
    
    # return the mass 
    return mass



# Printing Values for the Table

# MW Mass
MWHalo = ComponentMass("MW_000.txt",1)
MWDisk = ComponentMass("MW_000.txt",2)
MWBulge = ComponentMass("MW_000.txt",3)

print("MW Halo Mass", MWHalo, "(x 10^12)")
print("MW Disk Mass", MWDisk, "(x 10^12)" )
print("MW Bulge Mass", MWBulge, "(x 10^12)")
TotalMW = MWHalo + MWDisk + MWBulge
print("Total MW", TotalMW, "(x 10^12)")


# M31 Mass
M31Halo = ComponentMass("M31_000.txt",1)
M31Disk = ComponentMass("M31_000.txt",2)
M31Bulge = ComponentMass("M31_000.txt",3)

print("M31 Halo Mass", M31Halo, "(x 10^12)")
print("M31 Disk Mass", M31Disk, "(x 10^12)")
print("M31 Bulge Mass", M31Bulge, "(x 10^12)")
TotalM31 = M31Halo + M31Disk + M31Bulge
print("Total M31", TotalM31, "(x 10^12)")


# M33 Mass
M33Halo = ComponentMass("M33_000.txt",1)
M33Disk = ComponentMass("M33_000.txt",2)

print("M33 Halo Mass", M33Halo, "(x 10^12)")
print("M33 Disk Mass", M33Disk, "(x 10^12)")
TotalM33 = M33Halo + M33Disk
print("Total M33", TotalM33, "(x 10^12)")

print("Local Group Mass", TotalMW + TotalM31 + TotalM33, "(x 10^12)")




