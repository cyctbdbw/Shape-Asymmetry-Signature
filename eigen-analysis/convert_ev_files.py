# This script extracts the eigenvalues/eigenfunctions of the Shape-DNA output.
# See the references below for details of the Shape-DNA
# 1. Reuter, M. et al. Laplace–Beltrami spectra as ‘Shape-DNA’ of surfaces and solids. 
# Computer-Aided Design 38, 342-366 (2006).
# 2. Reuter, M. et al. Laplace-Beltrami Eigenvalues and Topological Features of 
# Eigenfunctions for Statistical Shape Analysis. Comput Aided Des 41, 739-755 (2009). 

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=FutureWarning) 

from argparse import ArgumentParser
import numpy as np


def convert_ev_files(evfile,out):
    #evfile='lh.sub1a.vtk-r0-d1-nbc.ev'
    area = ''
    volume = ''
    evals = []
    eigenValues = []
    ev_size1=[]
    ev_size2=[]
    # evecs = []
    with open(evfile, 'r') as inF:
        for line in inF:
            if 'Area:' in line:
                strlst = line.split()
                area = strlst[1]
            if 'Volume:' in line:
                strlst = line.split()
                volume = strlst[1]
            if 'Eigenvalues:' in line:
                evline = inF.next()
                evstr = ''
                while (evline is not None) and (not '}' in evline):
                    evstr = evstr+evline
                    evline = inF.next()
                evstr = evstr+evline
                evstr = evstr.translate(None,'{} \n')
                eigenValues = evstr.split(';')
                if abs(float(eigenValues[0])) < 10e-16 :
                    eigenValues[0] = "0"
                eigenValues_num=np.zeros([len(eigenValues),1])
                for i in range(0,len(eigenValues)):
                    eigenValues_num[i] = float(eigenValues[i])
                # ee=np.array([,1])
                # eigenValues.insert(0,volume)
                # eigenValues.insert(0,area)
            if 'Eigenvectors:' in line:
                size_line = inF.next()
                strlst=size_line.split()
                ev_size1=strlst[1]
                ev_size2=strlst[2]
                inF.next()
                evline = inF.next()
                evstr = ''
                ev_all = np.zeros([int(ev_size1),int(ev_size2)])
                while (evline is not None) and (not '}' in evline):
                    evstr = evstr+evline
                    evline = inF.next()
                evstr = evstr+evline
                evstr = evstr.translate(None,'{()} \n')
                evals = evstr.split(';')
                for i in range(0,int(ev_size2)):
                    ev_tmp=evals[i]
                    ev=ev_tmp.split(',')
                    for ni in range(0,int(ev_size1)):
                        value=float(ev[ni])
                        if abs(value) < (10e-16):
                            value=0
                        ev_all[ni,i]=float(ev[ni])
                        
    np.savetxt(out+'_evec.tsv', ev_all, delimiter="\t")
    np.savetxt(out+'_eval.tsv', eigenValues_num, delimiter="\t")

def main(raw_args=None):

    # Parse in inputs    
    parser = ArgumentParser(epilog="convert_ev_files.py -- A function to convert .ev files into .tsv files Kevin Aquino 2019 BMH")
    parser.add_argument("-ev", dest="evFile",
        help="Movement parameters", metavar="input.ev")
    parser.add_argument("-out", dest="out",
        help="The saved filed name for the Eigenvectors", metavar="out_ev")

    # import pdb;pdb.set_trace()
    # Here we are parsing the arguments
    args = parser.parse_args(raw_args)

    # Setting the arguments
    evFile = args.evFile
    out = args.out

    # Calculate FD here.
    convert_ev_files(evFile,out)

if __name__ == '__main__':
    main()
