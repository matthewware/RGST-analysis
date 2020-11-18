#!/opt/conda/bin/python
import matplotlib
matplotlib.use('Agg')

import datetime,time, pickle, os, sys, argparse

import numpy as np
from numpy.linalg import inv

import scipy as sp

def power_of_two(string):
    n = int(string)
    fl_lgn = np.floor(np.log(n)/np.log(2))
    if n != 2**fl_lgn:
        msg = "%r is not a power of two" % string
        raise argparse.ArgumentTypeError(msg)
    return n

parser = argparse.ArgumentParser(description='2 qubit GST analysis')
parser.add_argument('--maxlength', metavar='L', type=power_of_two, nargs=1, required=True,
                    help='Max germ power in experimental sequences (must be power of 2)')
parser.add_argument('--memory', metavar='M', type=float, nargs='?', default=10,
                    help='Memory available for analysis (in GB)')
parser.add_argument('--input', metavar='dir/data.gst', nargs=1, required=True,
                    help='Path to input data file')
parser.add_argument('--verbosity', metavar='V', type=int, nargs='?', required=True, default=1,
                    help='Verbosity level')
# parser.add_argument('--report', action='store_true',
#                     help='Generate pdf report')

args = parser.parse_args()

# This Logger objects ensures we keep a time-stamped log
# of the output of the pyGST analysis
import sys
class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        msg = message.replace("\n","\n"+datetime.datetime.now().strftime('%Y-%b-%d %H:%M:%S | '))
        self.terminal.write(msg)
        self.log.write(msg)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

import pygsti
from pygsti.objects import ConfidenceRegion
from pygsti.construction import std1Q_XYI

gs1Q = std1Q_XYI.gs_target
fiducials1Q = std1Q_XYI.fiducials
germs1Q = std1Q_XYI.germs

#maxLengths1Q = [0,1,2,4,8,16,32,64,128,256,512,1024]
maxLengths1Q = list(map(int,2**np.arange(0,int(np.log(args.maxlength)/np.log(2)+1))))

dsName = args.input[0]
pklName = "%s.pkl" % os.path.splitext(dsName)[0]
errPklName = "%s.errors.pkl" % os.path.splitext(dsName)[0]
sys.stdout = Logger(("%s."+datetime.datetime.now().strftime('%Y%b%d')+".log")% os.path.splitext(dsName)[0])

print()
print("Python version: ",sys.version)
print("Command line options")
print("  Memory limit (GB): ",args.memory)
print("  Max length: ",args.maxlength[0])
print("  Input file: ",args.input[0])
print()

start = time.time()
results      = pygsti.do_long_sequence_gst(dsName,
                                           gs1Q,
                                           fiducials1Q,
                                           fiducials1Q,
                                           germs1Q,
                                           maxLengths1Q,
                                           verbosity=3,
                                           advancedOptions={ 'memoryLimitInBytes' : args.memory*(1024)**3,
                                                             'depolarizeLGST' : 0.2} )
end = time.time()

# pickle dump before report generation, just in case
pickle.dump(results, open(pklName, "wb"))
print(("Written results to %s" % pklName))

hessian = pygsti.tools.logl_hessian(results.gatesets['final estimate'], results.dataset)
cr = ConfidenceRegion(gateset = results.gatesets['final estimate'], 
                      hessian = hessian,
                      confidenceLevel = 95,
                      hessianProjection = 'optimal gate CIs')

_gates = ['Gi','Gx','Gy']
def target_gate(i):
    return results.gatesets['target'][_gates[i]][:,:]

def mean_gate(i):
    return results.gatesets['final estimate'][_gates[i]][:,:]
    
def projected_hessian_for_gate(i):
    return (-cr.invRegionQuadcForm[(8+16*i):(8+16*(i+1)),(8+16*i):(8+16*(i+1))])

errors = {}
errors['proj_hessian'] = -cr.invRegionQuadcForm

errors['proj_gate_hessian'] = np.zeros((3,16,16))
errors['gate_cis'] = np.zeros((3,4,4))
errors['avg_fid'] = np.zeros(3)
errors['avg_fid_ci'] = np.zeros(3)

for idx in [0,1,2]:
    errors['proj_gate_hessian'][idx,:,:] = projected_hessian_for_gate(idx)
    cov = np.reshape(np.diag(errors['proj_gate_hessian'][idx,:,:]),(4,4))
    for row in range(4):
        for col in range(4):
            ci = sp.stats.norm.interval(0.95, 
                                        loc=results.gatesets['final estimate'][_gates[idx]][row,col], 
                                        scale=np.sqrt(cov[row,col]))
            errors['gate_cis'][idx,row,col] = (ci[1]-ci[0])/2.0
            
    fmean = (np.trace(mean_gate(idx) @ inv(target_gate(idx)))+2)/(4+2);
    fvar_p = target_gate(idx).flatten() @ projected_hessian_for_gate(2) @ target_gate(idx).flatten();
    ci_p = sp.stats.norm.interval(0.95, loc=fmean, scale=np.sqrt(fvar_p/(4.0+2.0)**2))
    errors['avg_fid'][idx] = fmean
    errors['avg_fid_ci'][idx] = (ci_p[1]-ci_p[0])/2.0
    
# pickle dump again, before report generation adds many fields to results object
pickle.dump(errors, open(errPklName, "wb"))
print(("Total time=%f hours" % ((end - start)/3600.0)))
