#!/usr/bin/env python
import numpy as np

def main():
    wmin=-6
    wmax=6
    nw=1000
    eta=0.02
    w = np.linspace(wmin,wmax,nw)

    fn = "g_coeffs.dat"
    f = open(fn,"r")
    nstep, norb = [int(x) for x in f.readline().split()]
    for ispin in range(1,3):
        for iorb in range(1,norb+1):
            print "Calculating (iorb,ispin)=(%d,%d)"%(iorb,ispin)

            fout = open("Aw.%d.%d.dat"%(iorb,ispin),"w")

            ev = float(f.readline().split()[2])

            Aw = np.zeros(nw, dtype=np.float64)

            ap = np.zeros(nstep,dtype=np.float64)
            bp = np.zeros(nstep,dtype=np.float64)
            an = np.zeros(nstep,dtype=np.float64)
            bn = np.zeros(nstep,dtype=np.float64)
            for istep in range(0,nstep):
                d1,d2,d3,d4 = [float(x) for x in f.readline().split()]
                ap[istep] = d1
                bp[istep] = d2
                an[istep] = d3
                bn[istep] = d4

            # continued fraction calculation
            i = 0
            for wp in w:
                z = wp+ev+eta*1j
                gr = continued_fraction_p(z, nstep, ap, bp)
                Aw[i] += -1/np.pi*gr.imag
                z = wp-ev+eta*1j
                gr = continued_fraction_m(z, nstep, ap, bp)
                Aw[i] += -1/np.pi*gr.imag
                i+=1
            i = 0
            for wp in w:
                fout.write("%e %e\n"%(wp,Aw[i]))
                i+=1
            fout.close()
    f.close()

def continued_fraction_p(z,nstep,ap,bp):

    f = bp[nstep-1]*bp[nstep-1]/(z-ap[nstep-1])
    for i in range(nstep-2,0,-1):
        f = bp[i]*bp[i]/(z-ap[i]-f)
    f = bp[0]*bp[0]/(z-ap[0]-f)
    return f

def continued_fraction_m(z,nstep,an,bn):

    f = bn[nstep-1]*bn[nstep-1]/(z+an[nstep-1])
    for i in range(nstep-2,0,-1):
        f = bn[i]*bn[i]/(z+an[i]-f)
    f = bn[0]*bn[0]/(z+an[0]-f)
    return f
if __name__ == "__main__":
    main()
