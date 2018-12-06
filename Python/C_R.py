import numpy as np
import math as m
import matplotlib.pyplot as plt
import time
import os

import readFunctions

def C_R(fname, inputFunctionHandle, I, cycleNo, noCycles, currentNo, noCurrents, LimUr, LimUc, pf1, pf2):
    n = eval(inputFunctionHandle)(fname, I)
    dataLen = n.shape[0]
    
    # Find voltage drop (discharge start)
    idx = -1
    startIdx = np.where(n[:,0] >= 0.5)[0][0]
    endIdx = np.where(n[:,0] >= 5)[0][0]
    base_du_dt = (n[endIdx,1] - n[startIdx,1]) / (n[endIdx,0] - n[startIdx,0])
    baseInterval = (n[endIdx,0] - n[startIdx,0]) / (endIdx - startIdx)
    found = False
    while not found:
        for i in range(startIdx,dataLen-1):
            if n[i + 1,0] - n[i,0] < baseInterval / 2:    # ignore points that are too close
                 continue
            du_dt = (n[i + 1,1] - n[i,1]) / (n[i + 1,0] - n[i,0])
            if du_dt < 5 * base_du_dt:     # discharge, so find 'more negative'
                idx = i
                Ur = n[idx,1]
                t0 = n[idx,0]             # discharge start index, voltage and time
                break
        if idx == -1:
            print("Unable to find voltage drop (cycle )" + str(cycleNo) + ")!")
            idx = 0
            Ur = n[idx,1]
            t0 = n[idx,0]             # discharge start index, voltage and time
            break
        # check for the false drop (voltage should decrease after the step for
        # at least 1s)
        endIdx_chk = np.where(n[:,0] >= t0 + 1)[0]
        if endIdx_chk.size == 0:    # if not found
            endIdx = dataLen - 1
        else:
            endIdx = endIdx_chk[0]

        for i in range(idx,endIdx):
            if n[i + 1,1] > n[i,1]:    # voltage raise
                break
        if i == endIdx - 1:                # voltage continuously decreasing
            found = True
        else:
            startIdx = idx + 1

    # Find end of discharge
    for i in range(idx,dataLen-3):
        if n[i+1,1] > n[i,1] and n[i+2,1] > n[i+1,1] and n[i+3,1] > n[i+2,1]:  # voltage raise for at least three samples
            break
    dataLen = i
    
    # Interval for calculating ESR
    if LimUr[0] == 1:    # 3rd order polynomial approximation between t0 and LimUr(2) is used to
                         # calculate ESR
        idxRStart = idx + 1
        idxREnd_chk = np.where(n[:,0] >= t0 + LimUr[1])[0]
        if idxREnd_chk.size == 0:    # if not found
            idxREnd = dataLen - 1
        else:
            idxREnd = idxREnd_chk[0]
        p = np.polyfit(n[idxRStart:idxREnd+1,0],n[idxRStart:idxREnd+1,1],3)   # 3rd order polynomial
    else:                   # linear approximation between voltages LimUr(1) and
                            # LimUr(2) is used to calculate ESR
        idxRStart_chk = np.where(n[:,1] < Ur * LimUr[0])[0]
        if idxRStart_chk.size == 0:  # if not found
            raise ValueError("Resistance calculation start voltage too low (cycle " + str(cycleNo) + ")!")
        else:
            idxRStart = idxRStart_chk[0]
        idxREnd_chk = np.where(n[:,1] < Ur * LimUr[1])[0]
        if idxREnd_chk.size == 0:  # if not found
            idxREnd = dataLen - 1
        else:
            idxREnd = idxREnd_chk[0]
        p = np.polyfit(n[idxRStart:idxREnd+1,0], n[idxRStart:idxREnd+1,1], 1)   # 1st order polynomial
        
    # Interval for calculating C
    if LimUc[0] == 1:           # Uc1 is voltage drop at the discharge start time
        idxCStart = np.where(n[:,1] < Ur)[0][0]   # find first data point having voltage < Ur

        if LimUc[1] == 0:        # Uc2 is voltage of last recorded discharge curve sample
            idxCEnd = dataLen - 1
        else:
            idxCEnd_chk = np.where(n[:,1] < Ur * LimUc[1])[0]
            if idxCEnd_chk.size == 0:  # if not found
                idxCEnd = dataLen - 1
            else:
                idxCEnd = idxCEnd_chk[0]
    else:
        idxCStart_chk = np.where(n[:,1] < Ur * LimUc[0])[0]
        if idxCStart_chk.size == 0:      # if not found
            raise ValueError("Capacitance calculation start voltage too low (cycle " + str(cycleNo) + ")!")
        else:
            idxCStart = idxCStart_chk[0]

        idxCEnd_chk = np.where(n[:,1] < Ur * LimUc[1])[0]
        if idxCEnd_chk.size == 0:  # if not found
            idxCEnd = dataLen - 1
        else:
            idxCEnd = idxCEnd_chk[0]

    UcStart = n[idxCStart,1]
    UcEnd = n[idxCEnd,1]


    # CAPACITANCE C = 2W / [(0.9ur)^2 - (0.7ur)^2], W - exchanged energy between UcStart and UcEnd
    # Energy
    W = 0
    for i in range(idxCStart+1, idxCEnd+1):
        W = W + abs(n[i,2]) * n[i,1] * (n[i,0]-n[i-1,0])    # dW = i * u * dt

    C = 2 * W / (UcStart*UcStart - UcEnd*UcEnd)

    # Self-discharge resistance
    U0 = n[0,1]
    Rd = - t0 / (C * m.log(Ur / U0))


    # ESR
    R = (Ur - np.polyval(p,t0)) / abs(n[idx+1,2])   # R = dU / I


    # Maximum power density
    Pd = 0.25 * Ur * Ur / R


    # Retrived energy
    idxStart = np.where(n[:,1] < Ur)[0][0]   # find first data point having voltage < Ur
    idxEnd = dataLen;

    E = 0
    for i in range(idxStart+1, idxEnd):
        E = E + abs(n[i,2]) * n[i,1] * (n[i,0]-n[i-1,0])    # dW = i * u * dt


    # NONLINEARITY
    # Ideal (linear) discharge curve
    t1x = n[idx,0]
    t1y = n[idx,1]

    t2x = n[dataLen-1,0]
    t2y = n[dataLen-1,1]

    # Nonlinearity (mean square deviation from ideal discharge curve
    k = (t2y - t1y) / (t2x - t1x)
    nl = 0
    for i in range(idx + 1, dataLen):
        x = n[i,0]
        y = t1y + k * (x - t1x)
        du = y - n[i,1]
        nl = nl + du * du
    nl = nl / (dataLen - idx)
    nl = m.sqrt(nl)


    # PLOTTING
    if cycleNo % pf1 == 0:      # plots are generated each 'pf1' cycle
        h = plt.figure(figsize = [10,6], dpi = 96)
        ax = h.add_subplot(1,1,1)
        ax.plot([t1x, t2x], [t1y, t2y], 'r--', linewidth = 1)
        ax.plot(n[:,0], n[:,1], 'b', linewidth = 2)
        plt.xlim([0, n[dataLen-1,0]])
        plt.ylim([n[dataLen-1,1], 1.05 * n[0,1]])

        # Markers
        ax.plot([0, n[idxCStart,0]], [n[idxCStart,1], n[idxCStart,1]],'m')
        my = n[idxCEnd,1]
        yl = plt.ylim()
        if my <= yl[0]:
            my = yl[0]+0.01
        ax.plot([0, n[idxCEnd,0]], [my, my], 'm')

        mEnd = n[idxREnd,0]+5  # end of ESR markers
        if mEnd > n[dataLen-1,0]:
            mEnd = n[dataLen-1,0]
        ax.plot([n[idxRStart,0], mEnd], [n[idxRStart,1], n[idxRStart,1]],'g')
        ax.plot([n[idxREnd,0], mEnd], [n[idxREnd,1], n[idxREnd,1]],'g')
    
        ax.plot([t0, t0], [np.polyval(p,t0), U0],':g', linewidth = 2)

        # discharge curve interpolation (linear or cubic)
        y = np.zeros(idxREnd-idx+1)
        for i in range(idx, idxREnd+1):
            y[i-idx] = np.polyval(p,n[i,0])
        ax.plot(n[idx:idxREnd+1,0], y,':c', linewidth = 2)
    
        ax.set_xlabel("t[s]")
        ax.set_ylabel("u[V]")

        plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

        txt = "Cycle: {cycleNo}\n\nC = {C:.3f} F\nESR = {R:.3f} \u03A9\nRd = {Rd:.0f} \u03A9"
        plt.text(0.8, 0.8, txt.format(cycleNo = cycleNo, C = C, R = R, Rd = Rd), bbox={'facecolor':'white', 'alpha':1, 'pad':10}, transform=ax.transAxes)

        plt.savefig(fname = os.path.splitext(fname)[0], dpi = 300)
        #plt.show(block=False)
        #plt.show()
        plt.close(h)
        print("Current: {I:.1f}\tCycle: {cn:d}".format(I = I*1000, cn = cycleNo))


    # Cumulative discharge plot
    if pf2 > 0 and cycleNo % pf2 == 0:
        noColors = noCycles / pf2
        col = (cycleNo / pf2 - 1) * 0.75 / noColors;  # R and G from 0 to 0.75

        if not hasattr(C_R, "cumulativeDischargePlot"):    # static variable, to hold figure handle
            C_R.cumulativeDischargePlot = plt.figure(figsize = [10,6], dpi = 96)
            
            ax1 = C_R.cumulativeDischargePlot.add_subplot(1,1,1)
            ax1.set_xlabel("t[s]")
            ax1.set_ylabel("u[V]")
            plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')
        else:
            plt.figure(num = C_R.cumulativeDischargePlot.number)

        plt.plot(n[0:dataLen,0], n[0:dataLen,1], color = [col, col, 1], linewidth = 2, label = str(cycleNo))
        plt.ylim([n[dataLen-1,1], 1.05 * n[0,1]])


    # Cumulative discharge plot by currents - for last cycle
    if not hasattr(C_R, "lastCurrentNo"):    # static variable, to determine when processing of next current started
        C_R.lastCurrentNo = -1
    
    if noCurrents > 1 and currentNo != C_R.lastCurrentNo and cycleNo == noCycles:
        C_R.lastCurrentNo = currentNo
        noColors = noCurrents
        col = currentNo * 0.75 / noColors;  # alter R and G from 0 to 0.75 for subsequent plots (shades of blue)
    
        if not hasattr(C_R, "dischargeIplot"):    # static variable, to hold figure handle
            C_R.dischargeIplot = plt.figure(figsize = [10,6], dpi = 96) # first plot, set up figure

            ax2 = C_R.dischargeIplot.add_subplot(1,1,1)
            ax2.set_xlabel("t[s]")
            ax2.set_ylabel("u[V]")
            plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')
        else:            # plot already exists
            plt.figure(num = C_R.dischargeIplot.number)

        plt.plot(n[0:dataLen,0],n[0:dataLen,1], color = [col, col, 1], linewidth = 2, label = str(1000*I) + " mA")
        plt.ylim([n[dataLen-1,1], 1.05 * n[0,1]])

    return C, R, Rd, nl, Pd, E

  