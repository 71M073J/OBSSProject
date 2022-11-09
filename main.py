import math

from scipy import ndimage
import wfdb
#import torch
import numpy as np
import numpy.lib.stride_tricks as stt
import matplotlib.pyplot as plt
import time
fs = 250
rec = wfdb.rdrecord("./long-term-st-database-1.0.0/s20011")
sig = rec.p_signal[:, 0][0:1250]#0:250 * 60 * 60]

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

#r_peaks = np.convolve(sig, filt)
r_peaks = np.convolve(sig, test, mode="same")
#plt.plot(r_peaks)
#plt.plot(r_peaks)
#plt.plot(r_peaks2)
#plt.show()
r_candidates = np.full(len(r_peaks), True)
left = (np.diff(r_peaks) >= 0)
right = (np.diff(r_peaks[::-1])[::-1] >= 0)
r_candidates[:-1] = r_candidates[:-1] & right
r_candidates[1:] = r_candidates[1:] & left
r_candidates = np.arange(len(r_peaks))[r_candidates]

second_order_diff = np.convolve(sig,np.array((-1, 2, -1)), mode="same")
c1 = 0.55
x2n = second_order_diff
sn = r_peaks * (sig + c1 * x2n)
plt.plot(sn, c="red")
plt.plot(x2n, c="green")
plt.plot(sig, c="magenta")
print(len(sn))
cands = []
#TODO OPTIMIZE THIS WIHT NUMPY
start_t = time.time()
maxl = len(sig) - 1
r = int(np.floor(fs * 0.2 + 0.5))
np.arange(len(sn))

#TODO mogoče loop do r da zračunamo vse offsete?
cands2 =[]
#np.abs(sn) < np.abs(sn - )
#np.where()
f = (np.abs(sn[r:-r]) < stt.sliding_window_view(np.abs(sn), 2 * r + 1).T).sum(axis=0)
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

asadasdasda = 0

end_t = time.time() + 1e-6
print(f"time elapsed for {len(sig)/fs} sec: ", end_t - start_t, f" speed {(len(sig)/fs)/(end_t - start_t)}s/s")
#TODO baseline extraction??? kaj je to, what do i do
#TODO MORE PREPROCESSING - we require accurate candidates
#print(r_candidates, cands)
#final = r_candidates
candidates = np.arange(len(sig))[cands]
final = set(r_candidates).intersection(set(cands))
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
plt.plot(r_peaks)
#plt.plot(second_order_diff)
#plt.plot(sig)
for i in final:
    ...
    plt.axvline(x=i, c="red")
plt.show()







quit()
#filter = np.cos(np.ones(20)/20 * 2 * np.pi)
#filter /= np.sum(filter)
#sig = np.convolve(sig, filter)
#fs = 250
print("pre-fft")
wave = np.fft.fft(sig)
print(("post-fft"))
#plt.plot(np.arange(len(wave)), np.real(wave))
#plt.plot(np.fft.fftfreq(wave.size, 1 / fs), wave.real)
# plt.plot(range(len(wave)), np.angle(wave), c="orange")
#plt.plot(range(len(wave)), np.imag(wave), c="red", alpha=0.3)
#plt.xlim(0)
print("post-plot")
#plt.show()
# print(vars(rec))
#print()
# wfdb.rdrecord()
# BRUH
# wtf is dduis
