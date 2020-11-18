#!/opt/conda/bin/python

import matplotlib
matplotlib.use('Agg')

import datetime, time, pickle, os, sys, argparse
import numpy as np

parser = argparse.ArgumentParser(description='Translation of GST PKL to npz data format')
parser.add_argument('--input', metavar='dir/analysis.pkl', nargs=1, required=True,
                    help='Path to analysis PKL file')

args = parser.parse_args()

import pygsti, pygsti.report

dsName = args.input[0]
pklName = "%s.pkl" % os.path.splitext(dsName)[0]
errPklName = "%s.errors.pkl" % os.path.splitext(dsName)[0]
pdfName = "%s.pdf" % os.path.splitext(dsName)[0]
npzName = "%s.npz" % os.path.splitext(dsName)[0]

print()
print("Python version: ",sys.version)
print("Command line options")
print("  Input file: ",args.input[0])
print()

result = pickle.load(open(pklName, 'rb'),encoding="latin1")
errors = pickle.load(open(errPklName, 'rb'),encoding="latin1")

#r.create_general_report_pdf(filename=pdfName,
#                            verbosity=args.verbosity)

final = np.zeros((3,4,4),dtype=complex)
target = np.zeros((3,4,4),dtype=complex)
pauli = np.zeros((3,4,4),dtype=complex)
cis = np.zeros((3,4,4),dtype=complex)
hessian = np.zeros((3,16,16),dtype=complex)
avgfid = np.zeros(3,dtype=complex)
avgfiderr = np.zeros(3,dtype=complex)

ptgs = result.gatesets['final estimate'].copy()

avgfid = errors['avg_fid']
avgfiderr = errors['avg_fid_ci']

for idx, g in enumerate(['Gi','Gx','Gy']):
    final[idx,:,:] = result.gatesets['final estimate'][g][:,:]
    target[idx,:,:] = result.gatesets['target'][g][:,:]
    hessian[idx,:,:] = errors['proj_gate_hessian'][idx,:,:]
    #cis[idx,:,:] = result.tables[u'bestGatesetGatesTable'].col(u'95% C.I. half-width')[idx]
    cis[idx,:,:] = errors['gate_cis'][idx,:,:]
    # pfinal is the idealized Pauli twirling of the final gateset,
    # assuming the target gates are Clifford group elements (multiplication is elemntwise)
    pauli[idx,:,:] = np.abs(target[idx,:,:])*final[idx,:,:]
    # Populate pauli-twirled gateset based on final estimate
    ptgs[g] = pauli[idx,:,:]

print(ptgs)

rows = len(result.tables['logLProgressTable']._rows)
logL = np.array([list(map(float,result.tables['logLProgressTable']._rows[r][0][0:8])) for r in range(rows)])

# Compute the likelihood of pauli twirled gateset (following pygsti/report/generation.py) for each subset of 'L'
# (similar to progress table, but without refitting at each step).
paulilogL = np.zeros((logL.shape[0],7))
for (idx, pair) in enumerate(zip(list(map(int, paulilogL[:,0])), result.gatestring_lists['iteration'])):
    L, gstrs = pair
    # print "%d, %d" % (L, len(gstrs))
    logL_upperbound = pygsti.tools.logl_max(result.dataset, gstrs)
    logl = pygsti.tools.logl(result.gatesets['final estimate'], result.dataset, gstrs)
    pt_logl = pygsti.tools.logl(ptgs, result.dataset, gstrs)
    if(logL_upperbound < logl):
        raise ValueError("LogL upper bound = %g but logl = %g!!" % (logL_upperbound, logl))
    if(logL_upperbound < pt_logl):
        raise ValueError("LogL upper bound = %g but Pauli-twirled logl = %g!!" % (logL_upperbound, pt_logl))
    Ns = len(gstrs)*(len(result.dataset.get_spam_labels())-1)
    Np = result.gatesets['final estimate'].num_nongauge_params()
    pt_Np = ptgs.num_nongauge_params()

    k = max(Ns-Np,0) #expected 2*(logL_ub-logl) mean
    pt_k = max(Ns-pt_Np,0) #expected 2*(logL_ub-logl) mean
    twoDeltaLogL = 2*(logL_upperbound - logl)
    pt_twoDeltaLogL = 2*(logL_upperbound - pt_logl)
    pt_Nsig = (pt_twoDeltaLogL-pt_k)/np.sqrt(2*pt_k)
    Nsig = (twoDeltaLogL-k)/np.sqrt(2*k)
    # put it in a table
    paulilogL[idx,0] = L
    paulilogL[idx,1] = twoDeltaLogL
    paulilogL[idx,2] = k
    paulilogL[idx,3] = Nsig
    paulilogL[idx,4] = pt_twoDeltaLogL
    paulilogL[idx,5] = pt_k
    paulilogL[idx,6] = pt_Nsig

# compute the log likelihood for some gateset in the context of some dataset
# [pygsti.logl( result.gatesets['final estimate'],
#               result.dataset,
#               result.gatestring_lists['iteration'][i])
#  for i in range(11)]

np.savez(npzName, 
         final=final, 
         cis=cis, 
         pauli=pauli, 
         target=target, 
         logL=logL, 
         paulilogL=paulilogL, 
         hessian=hessian,
         avgfid=avgfid,
         avgfiderr=avgfiderr)
