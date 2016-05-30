import numpy as np
from scipy import integrate


def main():

    # Set file locations
    totISRFfile   = "M31_212isrf_tot_ds_isrf.dat"
    oldISRFfile   = totISRFfile.replace('tot','old')
    yngISRFfile   = totISRFfile.replace('tot','yng')
    newISRFfile   = totISRFfile.replace('tot','new')
    
    # Get Labs values per component and per dust cell
    IDtot, x, y, z, Ltot = np.loadtxt(totISRFfile, usecols=(0,1,2,3,4,), unpack=True)
    IDold, Lold = np.loadtxt(oldISRFfile, usecols=(0,4,), unpack=True)
    IDyng, Lyng = np.loadtxt(yngISRFfile, usecols=(0,4,), unpack=True)
    IDnew, Lnew = np.loadtxt(newISRFfile, usecols=(0,4,), unpack=True)
    
    # write out
    writeLabsTot('Labs_tot.dat', IDtot, x, y, z, Ltot)
    writeLabs('Labs_old.dat', IDold, Lold)
    writeLabs('Labs_yng.dat', IDyng, Lyng)
    writeLabs('Labs_new.dat', IDnew, Lnew)

def writeLabsTot(file, ID, x, y, z, Labs):
    
    with open(file,'w') as f:
        f.write('# ID    x(pc)    y(pc)    z(pc)    Labs(W) \n')
        for id, xco, yco, zco, L in zip(ID, x, y, z, Labs):
            f.write(str(id)+'   '+str(xco)+'   '+str(yco)+'   '+str(zco)+'   '+str(L)+'\n')

def writeLabs(file, ID, Labs):
    
    with open(file,'w') as f:
        f.write('# ID    Labs(W) \n')
        for id, L in zip(ID,Labs):
            f.write(str(id)+'   '+str(L)+'\n')
    

if __name__ == '__main__':
    main()