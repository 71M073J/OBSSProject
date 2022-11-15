import collections
import math

from scipy import ndimage
import wfdb
#import torch
import numpy as np
import numpy.lib.stride_tricks as stt
import matplotlib.pyplot as plt
import time

fs = 250
filename = "./long-term-st-database-1.0.0/s20011"
anns = wfdb.rdann(filename, "atr")
true_beats = anns.sample
rec = wfdb.rdrecord(filename)
sig = rec.p_signal[:, 0][:]#0:250 * 60 * 60]

b1 = 0.025 * fs
b2 = 0.06 * fs
c = 2 * (b2 - b1) / (2 * b1 + 1)
filt_size = 51
filt = np.zeros(filt_size)
print(b1, b2, c)
for i in range(-len(filt)//2,len(filt)//2,1):
    if -b1 < i < b1:
        filt[i + len(filt)//2] = c
    elif b1 < np.abs(i) <= b2:
        filt[i + len(filt)//2] = -1
#plt.plot(range(-len(filt)//2,len(filt)//2,1),filt)
impulse = np.zeros(filt_size)
impulse[len(filt)//2] = 1
test = 6 * fs * ndimage.gaussian_laplace(-impulse, b1) # TESTNI FILTER ZA LoG
#plt.plot(range(-len(filt)//2,len(filt)//2,1),test)
#plt.show()
#test = filt#TODO a bomo LoG uporablal al ne
#r_peaks = np.convolve(sig, filt)
r_peaks = np.convolve(sig, test, mode="same")
#plt.plot(r_peaks)
#plt.plot(r_peaks)
#plt.plot(r_peaks2)
#plt.show()
#r_candidates = np.full(len(r_peaks), True)
#left = (np.diff(r_peaks) >= 0)
#right = (np.diff(r_peaks[::-1])[::-1] >= 0)
#r_candidates[:-1] = r_candidates[:-1] & right
#r_candidates[1:] = r_candidates[1:] & left
#r_candidates = np.arange(len(r_peaks))[r_candidates]

second_order_diff = np.convolve(sig,np.array((-1, 2, -1)), mode="same")
c1 = 0.55
x2n = second_order_diff
sn = r_peaks * (sig + c1 * x2n)
#plt.plot(sn, c="red")
#plt.plot(x2n, c="green")
#plt.plot(sig, c="magenta")
print(len(sn))
cands = []

maxl = len(sig) - 1
r = int(np.floor(fs * 0.2 + 0.5))
np.arange(len(sn))

#np.abs(sn) < np.abs(sn - )
#np.where()
abssn = np.abs(sn)
f = (np.abs(sn[r:-r]) < stt.sliding_window_view(abssn, 2 * r + 1).T).sum(axis=0)
cands = np.arange(r, len(f) + r)[f == 0]
#f = f.sum(axis=0)
#css = [[] for i in sn]
#for k, i in enumerate(sn):
#    fail = False
#    if not (np.abs(sn[k]) < np.abs(sn[np.clip(int(k - r), a_min=0, a_max=maxl):np.clip(int(k + r + 1), a_min=0, a_max=maxl)])).any():
#        cands2.append(k)
#    css[k].append(
#        np.abs(sn[k]) < np.abs(sn[
#                               np.clip(int(k - r), a_min=0, a_max=maxl):
#                               np.clip(int(k + r), a_min=0, a_max=maxl)]))

start_t = time.time()
T = 0.6
# np.argpartition([...])#maybe celo hitreje da gremo čez sn namesto čez kandidate
#s5 = np.argpartition(stt.sliding_window_view(np.abs(sn), 10 * fs), -5,axis=1)[:,-5:]
#NVM this array is 384 GiB ()
w1 = np.zeros(len(cands))
w2 = np.zeros(len(cands))
taus = np.arange(1,6).astype(np.float32)[::-1]
taus /= np.sum(taus)
for ind, i in enumerate(cands):#navadn loop bo še najhitrej
    small = False
    curind = i
    idx = ind
    max5 = np.zeros(5)
    minind = 0
    while idx > -1 and (curind:=cands[idx]) + 10 * fs > i:
        #curind = cands[idx]
        curr = abssn[curind]
        if ind - idx < 5:
            minind = np.argmin(max5)
        else:
            minind = np.nonzero(max5)[0][np.argmin(max5[np.nonzero(max5)])]
        if curr > max5[minind]:
            max5[minind] = curr
        #get max 5
        idx -= 1
    w1[ind] = max5[np.nonzero(max5)[0][np.argmin(max5[np.nonzero(max5)])]]

    ms = np.zeros(7)
    mind = np.clip(ind-6, a_min=0, a_max=ind)
    ms[1:ind - mind + 1] = cands[mind:ind][::-1]
    ms[0] = i
    difms = np.diff(ms[::-1])[::-1]
    if difms[1:].nonzero()[0].any() and difms[0] < np.mean(difms[1:][difms[1:].nonzero()]) * 0.7:
        difms[0:5] = difms[1:6]
        small = True
    num_nonzero = difms[:-1].nonzero()[0].shape[0]
    ie = np.sum((taus[:num_nonzero] / np.sum(taus[:num_nonzero])) * difms[:num_nonzero])
    beta1 = 1
    beta2 = 1
    if small:
        w2[ind] = beta1 + beta2 * np.abs(((ms[0] - ms[2])/2) / ie - np.floor(0.5 + (ms[0] - ms[1]) / ie))
    else:
        w2[ind] = beta1 + beta2 * np.abs((ms[0] - ms[1])/ie - np.floor(0.5 + (ms[0] - ms[1])/ie))


w1 *= T

adaptive_threshold = w1 * w2
adapted_cands = cands[adaptive_threshold < np.abs(sn[cands])]
end_t = time.time() + 1e-6

print(adaptive_threshold, np.abs(sn[cands]))
print(len(cands), "dolžina kandidatov")
print(f"time elapsed for {len(sig)/fs} sec: ", end_t - start_t, f" speed {(len(sig)/fs)/(end_t - start_t)}s/s")
print(len(true_beats), len(adapted_cands), "KEKW")
#plt.plot(cands, adaptive_threshold, c="red")
#plt.plot(sig)
#plt.show()
#TODO baseline extraction??? kaj je to, what do i do
#TODO MORE PREPROCESSING - we require accurate candidates
#print(r_candidates, cands)
#final = r_candidates
candidates = np.arange(len(sig))[cands]
#final = set(r_candidates).intersection(set(cands))
final = cands
#final = [r_candidates[0]]
#for i in range(1, len(r_candidates)):
#    if r_candidates[i] > final[-1] + fs/5:
#        final.append(r_candidates[i])

#print(sig)
#print(np.diff(sig))
#print(np.diff(np.diff(sig)))
#final = cands + r + 1
#print(final)
plt.figure(figsize=(15,5))
plt.plot(r_peaks)
#plt.plot(second_order_diff)
#plt.plot(sig)
for i in final:
    plt.axvline(x=i, c="red")
for i in adapted_cands:
    print(i)
    plt.axvline(x=i, c="blue")
plt.show()
