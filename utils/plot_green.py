#!/usr/bin/env python
f=open("green.save")
na, nspin, norb, nw = [int(x) for x in f.readline().split()]

nwcut = 100

for ia in range(1,na+1):
    for ispin in range(1,nspin+1):
        for iorb in range(1,norb+1):
            f.readline()
            fn = "green.%d.%d.%d.dat"%(ia,ispin,iorb)
            print fn
            fout=open(fn,"w")
            for iw in range(0,nw):
                l = f.readline()
                if iw<nwcut:
                    fout.write(l)
            fout.close()


