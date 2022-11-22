import collections
import math

import matplotlib.animation as animation
from scipy import ndimage
import wfdb
# import torch
import numpy as np
import numpy.lib.stride_tricks as stt
import matplotlib.pyplot as plt
import time
from scipy.signal import firwin, firwin2, freqz
import json


def detect(sig, true_beats, T=0.7, beta1=1.0, beta2=1.0, fs=250, search=False, name="default",
           maxbeta=2.0, minT=0.5, maxT=1.0, signum=0, annotate=False):
    #print(len(sig))
    #TODO BASELINE EXTRACTION
    # HIGH_PASS FILTER
    # FIR FILTER
    # cutoff okoli 0.8 Hz
    t0 = time.time()
    cutoff_Hz = 0.8
    nyq_freq = fs / 2
    fir_filter_highpass08 = firwin(1001, cutoff_Hz/nyq_freq, pass_zero="highpass")
    #plt.plot(sig)
    sig = np.convolve(sig, fir_filter_highpass08, mode="same")
    sig -= np.mean(sig)
    #plt.show()
    print(f"Baseline Extraction duration:{time.time() - t0} seconds")
    #quit()
    t1 = time.time()
    b1 = 0.025 * fs
    b2 = 0.06 * fs
    c = 2 * (b2 - b1) / (2 * b1 + 1)
    filt_size = 51
    filt = np.zeros(filt_size)
    #print(b1, b2, c)
    for i in range(-len(filt) // 2, len(filt) // 2, 1):
        if -b1 < i < b1:
            filt[i + len(filt) // 2] = c
        elif b1 < np.abs(i) <= b2:
            filt[i + len(filt) // 2] = -1
    # plt.plot(range(-len(filt)//2,len(filt)//2,1),filt)
    impulse = np.zeros(filt_size)
    impulse[len(filt) // 2] = 1
    test = 6 * fs * ndimage.gaussian_laplace(-impulse, b1)  # TESTNI FILTER ZA LoG
    # plt.plot(range(-len(filt)//2,len(filt)//2,1),test)
    # plt.show()
    # test = filt#TODO a bomo LoG uporablal al ne
    # r_peaks = np.convolve(sig, filt)
    r_peaks = np.convolve(sig, test, mode="same")
    print(f"LoG filtering duration: {time.time() - t1} seconds.")
    # plt.plot(r_peaks)
    # plt.plot(r_peaks)
    # plt.plot(r_peaks2)
    # plt.show()
    # r_candidates = np.full(len(r_peaks), True)
    # left = (np.diff(r_peaks) >= 0)
    # right = (np.diff(r_peaks[::-1])[::-1] >= 0)
    # r_candidates[:-1] = r_candidates[:-1] & right
    # r_candidates[1:] = r_candidates[1:] & left
    # r_candidates = np.arange(len(r_peaks))[r_candidates]

    second_order_diff = np.convolve(sig, np.array((-1, 2, -1)), mode="same")
    c1 = 0.55
    x2n = second_order_diff
    sn = r_peaks * (sig + c1 * x2n)
    # plt.plot(sn, c="red")
    # plt.plot(x2n, c="green")
    # plt.plot(sig, c="magenta")
    cands = []

    start_t = time.time()
    maxl = len(sig) - 1
    r = int(np.floor(fs * 0.2 + 0.5))
    np.arange(len(sn))


    abssn = np.abs(sn)
    f = (np.abs(sn[r:-r]) < stt.sliding_window_view(abssn, 2 * r + 1).T).sum(axis=0)
    cands = np.arange(r, len(f) + r)[f == 0]
    print(f"Candidate filtering 1 duration: {time.time() - start_t} seconds.")
    # f = f.sum(axis=0)
    # css = [[] for i in sn]
    # for k, i in enumerate(sn):
    #    fail = False
    #    if not (np.abs(sn[k]) < np.abs(sn[np.clip(int(k - r), a_min=0, a_max=maxl):np.clip(int(k + r + 1), a_min=0, a_max=maxl)])).any():
    #        cands2.append(k)
    #    css[k].append(
    #        np.abs(sn[k]) < np.abs(sn[
    #                               np.clip(int(k - r), a_min=0, a_max=maxl):
    #                               np.clip(int(k + r), a_min=0, a_max=maxl)]))
    # np.argpartition([...])#maybe celo hitreje da gremo čez sn namesto čez kandidate
    # s5 = np.argpartition(stt.sliding_window_view(np.abs(sn), 10 * fs), -5,axis=1)[:,-5:]
    # NVM this array is 384 GiB ()
    t3 = time.time()
    w1 = np.zeros(len(cands))
    w2 = np.zeros(len(cands))
    taus = np.arange(1, 6).astype(np.float32)[::-1]
    taus /= np.sum(taus)
    for ind, i in enumerate(cands):  # navadn loop bo še najhitrej
        small = False
        #curind = i
        idx = ind
        max5 = np.zeros(5)
        #minind = 0
        while idx > -1 and (curind := cands[idx]) + 10 * fs > i:
            # curind = cands[idx]
            curr = abssn[curind]
            if ind - idx < 5:
                minind = np.argmin(max5)
            else:
                minind = np.nonzero(max5)[0][np.argmin(max5[np.nonzero(max5)])]
            if curr > max5[minind]:
                max5[minind] = curr
            # get max 5
            idx -= 1
        w1[ind] = max5[np.nonzero(max5)[0][np.argmin(max5[np.nonzero(max5)])]]

        ms = np.zeros(7)
        mind = np.clip(ind - 6, a_min=0, a_max=ind)
        ms[1:ind - mind + 1] = cands[mind:ind][::-1]
        ms[0] = i
        difms = np.diff(ms[::-1])[::-1]
        if difms[1:].nonzero()[0].any() and difms[0] < np.mean(difms[1:][difms[1:].nonzero()]) * 0.7:
            difms[0:5] = difms[1:6]
            small = True
        num_nonzero = difms[:-1].nonzero()[0].shape[0]
        ie = np.sum((taus[:num_nonzero] / np.sum(taus[:num_nonzero])) * difms[:num_nonzero])
        if small:
            w2[ind] = np.abs(((ms[0] - ms[2]) / 2) / ie - np.floor(0.5 + (ms[0] - ms[1]) / ie))
        else:
            w2[ind] = np.abs((ms[0] - ms[1]) / ie - np.floor(0.5 + (ms[0] - ms[1]) / ie))
    print(f"Adaptive threshold filtering duration: {time.time() - t3} seconds.")
    if search:
        w1save = np.copy(w1)
        w2save = np.copy(w2)

        maxscore = len(true_beats)
        true = set(true_beats)
        maxs = 0
        maxvals = 0
        a = np.arange(0.0, 0.4, 0.05)
        b = np.arange(0., maxbeta, 0.25)
        c = np.arange(0., maxbeta, 0.25)
        ims = np.zeros((len(a), len(b), len(c)))
        t_s = time.time()
        maxvalslist = []
        for i1, T in enumerate(a):
            for i2, beta1 in enumerate(b):
                for i3, beta2 in enumerate(c):
                    w1 = w1save * T
                    w2 = beta1 + beta2 * w2save
                    adaptive_threshold = w1 * w2
                    adapted_cands = cands[adaptive_threshold < np.abs(sn[cands])]
                    # scoring:
                    tp = 0
                    fp = 0
                    fn = 0
                    cand = set(adapted_cands)
                    for can in cand:
                        gotOne = False
                        for offset in range(-5,6,1):
                            if can + offset in true:
                                gotOne = True
                                break
                        if gotOne:
                            tp += 1
                        else:
                            fp += 1
                    for tr in true:
                        gotOne = False
                        for offset in range(-5, 6, 1):
                            if tr + offset in cand:
                                gotOne = True
                                break
                        if not gotOne:
                            fn += 1
                    print(f"True Positive: {tp/maxscore}, \nFalse Positive:"
                          f" {fp/maxscore} \nFalse Negative: {fn/maxscore}")
                    print("params", T, beta1, beta2,"maxs", maxs, len(adapted_cands), maxscore)
                    #if tp/(fp + 1) > maxs:
                    #    maxs = tp/(fp + 1)
                    #    maxvals = (T, beta1, beta2)
                    fs = tp / (tp + 0.5 * (fp + fn))
                    ims[i1, i2, i3] = fs
                    print("F-score:", fs)
                    maxvalslist.append((fs,(T, beta1, beta2)))
                    if fs > maxs:
                        maxs = fs
                        maxvals = (T, beta1, beta2)
        maxvalslist = [x[1] for x in maxvalslist if x[0] == maxs]
        print(maxs, maxvals)
        print(f"Search duration: {time.time() - t_s} seconds.")
        with open("heartbeat_search.json", "rb") as f:
            dicc = json.load(f)
            dicc[name + str(signum)] = maxvals
            dicc[name + str(signum) + "f_score"] = maxs
            dicc[name + str(signum) + "maxvals"] = maxvalslist
        with open("heartbeat_search.json", "w") as f:
            json.dump(dicc, f)
        fig, ax = plt.subplots()
        anims = []
        for h in ims:
            im = ax.imshow(h, animated=True, vmin=ims.min(), vmax=ims.max())
            anims.append([im])
        ani = animation.ArtistAnimation(fig, anims, interval=500, blit=True,
                                        repeat_delay=1000)
        ani.save(f"./{name}{signum}.mp4")
        #plt.show()
    else:
        w1 *= T
        w2 = beta1 + beta2 * w2
        adaptive_threshold = w1 * w2
        adapted_cands = cands[adaptive_threshold < np.abs(sn[cands])]
        end_t = time.time() + 1e-6
        cand = set(adapted_cands)
        tp = 0
        fp = 0
        fn = 0
        maxscore = len(true_beats)
        true = set(true_beats)
        maxs = 0
        for can in cand:
            gotOne = False
            for offset in range(-5, 6, 1):
                if can + offset in true:
                    gotOne = True
                    break
            if gotOne:
                tp += 1
            else:
                fp += 1
        for tr in true:
            gotOne = False
            for offset in range(-5, 6, 1):
                if tr + offset in cand:
                    gotOne = True
                    break
            if not gotOne:
                fn += 1
        print(f"True Positive: {tp / maxscore}, \nFalse Positive:"
              f" {fp / maxscore} \nFalse Negative: {fn / maxscore}")
        print("params", T, beta1, beta2, "maxs", maxs, len(adapted_cands), maxscore)
        fs = tp / (tp + 0.5 * (fp + fn))
        print("F-score:", fs)
        print(f"time elapsed for {len(sig) / fs} sec: ", end_t - t0,
              f" speed {(len(sig) / fs) / (end_t - start_t)}s/s")
        sensitivity = tp / (tp + fn)
        pos_pred = tp / (tp + fp)
        if annotate:
            with open("det/" + name + ".det", "w") as f:
                for i in adapted_cands:
                    f.write(i)

    #print(adaptive_threshold, np.abs(sn[cands]))
    #print(len(cands), "dolžina kandidatov,")
    #print(len(true_beats), len(adapted_cands), "resničnih QRS in pa končnih kandidatov")


if __name__ == "__main__":
    fs = 250
    for h in range(1,10):
        filename = f"./long-term-st-database-1.0.0/s200{h}1"
        anns = wfdb.rdann(filename, "atr")
        true_beats = anns.sample[anns.sample < 1250]
        true_beats = anns.sample
        rec = wfdb.rdrecord(filename)

        for i in range(1):
            sig = rec.p_signal[:, i][:]  # 0:250 * 60 * 60]
            #detect(sig, true_beats, search=True, name=filename.split("/")[-1], maxbeta=2.6, signum=i)

            detect(sig, true_beats, T=0.15, beta1=1., beta2=1.)
    #detect(sig, true_beats, T=0.68, beta1=0.6, beta2=1.8)
    #print(f"T: {T}, b1: {b1}, b2: {b2}, max score: {maxs}")


quit()
cands = ...
r_peaks = ...
adapted_cands = ...
# TODO baseline extraction??? kaj je to, what do i do
# TODO MORE PREPROCESSING - we require accurate candidates
candidates = np.arange(len(sig))[cands]
final = cands
plt.figure(figsize=(15, 5))
plt.plot(r_peaks)
# plt.plot(second_order_diff)
# plt.plot(sig)
for i in final:
    plt.axvline(x=i, c="red")
for i in adapted_cands:
    plt.axvline(x=i, c="blue")
plt.show()
