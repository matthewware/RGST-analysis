#!/opt/conda/bin/python

import datetime,time, pickle, os, sys, argparse

import numpy as np

def power_of_two(string):
    n = int(string)
    fl_lgn = np.floor(np.log(n)/np.log(2))
    if n != 2**fl_lgn:
        msg = "%r is not a power of two" % string
        raise argparse.ArgumentTypeError(msg)
    return n

parser = argparse.ArgumentParser(description='1 qubit GST experiment generation')
parser.add_argument('--maxlength', metavar='L', type=power_of_two, nargs=1, required=True,
                    help='Max germ power in experimental sequences (must be power of 2)')

args = parser.parse_args()

import pygsti
from pygsti.construction import std1Q_XYI

gs1Q = std1Q_XYI.gs_target
fiducials1Q = std1Q_XYI.fiducials
germs1Q = std1Q_XYI.germs
#effects1Q = std1Q_XYI.effectStrs

#maxLengths1Q = [0,1,2,4,8,16,32,64,128,256,512,1024]
maxLengths1Q = list(map(int,2**np.arange(0,int(np.log(args.maxlength)/np.log(2)+1))))

listOfExperiments = pygsti.construction.make_lsgst_experiment_list( gateLabels = list(gs1Q.gates.keys()),
                                                                    prepStrs = fiducials1Q,
                                                                    effectStrs = fiducials1Q,
                                                                    germList = germs1Q,
                                                                    maxLengthList = maxLengths1Q)


pygsti.io.write_gatestring_list('gst-1q-1024.txt', listOfExperiments)
