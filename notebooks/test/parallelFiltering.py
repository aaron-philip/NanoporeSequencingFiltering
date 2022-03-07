from multiprocessing import Queue, Process
from filterPy import TimesteppedFilters, ButterworthFilters
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy import signal
from scipy.optimize import curve_fit
import math
from typing import Optional
import os




def queueFilters(obj, truth, ftype, q, order=None):
    if order is not None:
        q.put(obj.timestepped.optSolve(order, truth, ftype=ftype))
    else:
        q.put(obj.timestepped.optSolve(truth, ftype=ftype))



def parallelizeCalcs(): # spawns child processes
    q = Queue()
    processes = []
    rets = []
    p1 = Process(target=queueFilters, args=(timestepped, subTruth[sample], 'lowpass', q))
    processes.append(p1)
    p1.start()
    p2 = Process(target=queueFilters, args=(timestepped, subTruth[sample], 'highpass', q))
    processes.append(p2)
    p2.start()
    p3 = Process(target=queueFilters, args=(timestepped, subTruth[sample], 'bandpass', q))
    processes.append(p3)
    p3.start()
    p4 = Process(target=queueFilters, args=(butterworth, subTruth[sample], 'lowpass', q, 10))
    processes.append(p4) 
    p4.start()

    for p in processes:
        ret = q.get() # will block
        rets.append(ret)
    for p in processes:
        p.join()
    return rets

if __name__== "__main__":
    print(parallelizeCalcs())
