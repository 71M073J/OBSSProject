import collections
import math
import os

import matplotlib.animation as animation
from scipy import ndimage
import wfdb
# import torch
import numpy as np
import numpy.lib.stride_tricks as stt
import matplotlib.pyplot as plt
import time
from scipy.signal import firwin, firwin2, freqz
from scipy.optimize import minimize
import json


def detect(sig, true_beats, signum, T=0.7, beta1=1.0, beta2=1.0, fs=250, search=False, name="default",
           annotate=False, verbose=True, plot=False, Tlim=(5.0, 5.0), beta1lim=(0.0, 1.0), beta2lim=(0.0, 1.0)):
    # print(len(sig))
    # HIGH_PASS FILTER
    # FIR FILTER
    # cutoff okoli 0.8 Hz
    t0 = time.time()
    cutoff_Hz = 0.8
    nyq_freq = fs / 2
    fir_filter_highpass08 = firwin(1001, cutoff_Hz / nyq_freq, pass_zero="highpass")
    if plot and len(sig) > fs * 60 * 5:
        print("Will not plot over 5 min of data; setting plot to false.")
        # plot = False
    if plot:
        plt.plot(sig, c="purple")
    sig = np.convolve(sig, fir_filter_highpass08, mode="same")
    sig -= np.mean(sig)
    sig /= np.abs(sig).max()
    if plot:
        plt.plot(sig, c="red")
    if verbose:
        print(f"Baseline Extraction duration:{time.time() - t0} seconds")
    # quit()
    t1 = time.time()
    b1 = 0.025 * fs
    b2 = 0.06 * fs
    c = 2 * (b2 - b1) / (2 * b1 + 1)
    filt_size = 51
    filt = np.zeros(filt_size)
    # print(b1, b2, c)
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
    if verbose:
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
    if verbose:
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
    maxs = 0
    for ind, i in enumerate(cands):  # navadn loop bo še najhitrej
        small = False
        # curind = i
        idx = ind
        max5 = np.zeros(5)
        # minind = 0
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
    if verbose:
        print(f"Adaptive threshold filtering duration: {time.time() - t3} seconds.")
    if search:
        w1save = np.copy(w1)
        w2save = np.copy(w2)

        maxscore = len(true_beats)
        true = set(true_beats)
        maxs = 0

        def eval_params(x):
            T, beta1, beta2 = x
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

            # if tp/(fp + 1) > maxs:
            #    maxs = tp/(fp + 1)
            #    maxvals = (T, beta1, beta2)
            fs = tp / (tp + 0.5 * (fp + fn))
            # if verbose:
            sensitivity = tp / (tp + fn + 1e-12)
            pos_pred = tp / (tp + fp + 1e-12)
            score = (sensitivity + pos_pred)
            print("Score:", score / 2)
            return -fs
            # ims[i1, i2, i3] = fs
            # print("F-score:", fs)
            # maxvalslist.append((fs, (T, beta1, beta2)))

        res = minimize(eval_params, np.array((T, beta1, beta2)), method="Nelder-Mead",
                       bounds=[Tlim, beta1lim, beta2lim])
        print(res.x)
        with open("min_params.txt", "a") as fa:
            fa.write(str(res.x) + ", F-score: " + str(res.fun))
    else:
        print(np.mean(sn[cands]))
        w1 *= T
        w2 = beta1 + beta2 * w2
        adaptive_threshold = w1 * w2
        adapted_cands = cands[adaptive_threshold < np.abs(sn[cands])]
        cand = set(adapted_cands)
        tp = 0
        fp = 0
        fn = 0
        maxscore = len(true_beats)
        true = set(true_beats)
        maxs = 0
        for can in cand:
            gotOne = False
            for offset in range(-int(fs * 0.15 / 2), int(fs * 0.15 / 2) + 1, 1):
                if can + offset in true:
                    gotOne = True
                    break
            if gotOne:
                tp += 1
            else:
                fp += 1
        for tr in true:
            gotOne = False
            for offset in range(-int(fs * 0.15 / 2), int(fs * 0.15 / 2) + 1, 1):
                if tr + offset in cand:
                    gotOne = True
                    break
            if not gotOne:
                fn += 1

        end_t = time.time() + 1e-6
        if verbose:
            print(f"True Positive: {tp / maxscore}, \nFalse Positive:"
                  f" {fp / maxscore} \nFalse Negative: {fn / maxscore}")
            print("params", T, beta1, beta2, "maxs", maxs, len(adapted_cands), maxscore)
        f_s = tp / (tp + 0.5 * (fp + fn))
        print("F-score:", f_s)
        print(f"time elapsed for {len(sig) / fs} sec: ", end_t - t0,
              f"sec, speed {(len(sig) / fs) / (end_t - start_t)} s/s")
        sensitivity = tp / (tp + fn + 1e-12)
        pos_pred = tp / (tp + fp + 1e-12)
        print(f"Sensitivity: {sensitivity}\n+Predictivity: {pos_pred}")
        if plot:
            for x in adapted_cands:
                plt.axvline(x, alpha=0.5)
            plt.show()
        if annotate:
            if not os.path.exists("./det"):
                os.mkdir("./det")
            with open(f"det/{name}_signal{signum}.det", "w") as f:
                data = "\n".join([f"0:00:00.00 {x} N 0 0 0" for x in adapted_cands])
                f.write(data)
        return sensitivity, pos_pred


if __name__ == "__main__":
    # TODO use minimiser instead of search
    fs = 250
    prefix = "E:/Users/timotej999/Downloads/long-term-st-database-1.0.0/"
    results = []
    db = f"{prefix}./long-term-st-database-1.0.0/"
    done = set()
    for p in os.listdir(db):
        if os.path.isdir(p):
            continue
        if os.path.exists(db + p.split("/")[-1].split(".")[0] + ".hea"):
            filename = db + p.split("/")[-1].split(".")[0]
        else:
            continue
        if filename.split("/")[-1] in done:
            continue
        # for h in range(1, 10):
        #    filename = f"{prefix}./long-term-st-database-1.0.0/s200{h}1"
        # filename = f"{prefix}./long-term-st-database-1.0.0/s20501"
        anns = wfdb.rdann(filename, "atr")
        # true_beats = anns.sample[anns.sample < 1250]
        true_beats = anns.sample
        rec = wfdb.rdrecord(filename)

        for i in range(2):
            sig = rec.p_signal[:, i][:]  # [34*60*250:35*60*250]  # 0:250 * 60 * 60]

            # plt.plot(sig[34*60*250:35*60*250])
            # plt.show()
            search = False
            sens, pluspred = detect(sig, true_beats, search=search, name=filename.split("/")[-1], T=0.3, beta1=0.9,
                                    beta2=1.3, signum=i, verbose=True, Tlim=(0.0, 2.0), beta1lim=(0.0, 2.),
                                    beta2lim=(0.0, 2.), plot=False, annotate=True)
            if i == 0:
                results.append([filename.split("/")[-1], sens, pluspred])
            else:
                results[-1] += [sens, pluspred]
        done.add(filename.split("/")[-1])
    with open("./det/all.txt", "w") as f:
        print("FILENAME, s1 SENS, s1 +PRED, s2 SENS, s2 +PRED", file=f)
        for i in results:
            print(f"{','.join([str(x) for x in i])}", file=f)
    # detect(sig, true_beats, T=0.68, beta1=0.6, beta2=1.8)
    # print(f"T: {T}, b1: {b1}, b2: {b2}, max score: {maxs}")
