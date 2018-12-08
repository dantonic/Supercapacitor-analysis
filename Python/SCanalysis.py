import glob     # UNIX style pathname pattern expansion
import os
import numpy as np
import matplotlib.pyplot as plt

import C_R
import annotate

def SCanalysis(Ur, I, m, S, text1, text2, text3, text4, text5, files, dataFunct, LimUr, LimUc, pf1,pf2):
    """
    SCanalysis Calculate supercapacitor capacitance and equivalent serial resistance

    SCanalysis(I,Ur,m,S,'text1','text2','text3','text4','text5',files,[Ur1 Ur2],[Uc1 Uc2],pf1,pf2)

    Parameters:
       Ur[V] - Nominal voltage (to which SC is charged)
       I [mA] - Discharge current
           a)	Scalar (e.g.21.5) Analysis for single current is performed
           b)	Vector (e.g.[10 17.8 31.6 56.2 100]) Performs multi cycle - multi current analysis.
       m[mg] - Mass of active material. Should be zero if not used.
       S[cm2] - Area of one electrode. Should be zero if not used.
       text1..text5 - Five comment lines that will be printed at plots
       files - File name filter.Files that match the filter are sorted according to the file name.
               Each file should contain data for one discharge cycle.
               Care should be taken to ensure that sorting order corresponds to the cycles order,
               e.g.name files as '17_04_25_Data_001.txt' ... '17_04_25_Data_100.txt'.In this case
               'fnames' field should be set to *Data*.txt,where '*' matches any number of characters.
               If vector of discharge currents is defined,files field should be defined as cell array,
               having same number of elements as the I vector (and optionally one additional element),
               where each element defines file name filter for corresponding current measurements,
                   e.g.: {‘I_15mA*.txt’ ‘I_22mA*.txt’}
               If additional file is defined, it should contain repeated measurement for the first current.
               Data from that file is shown with red cross at plots.
       dataFunct - Handle of function for reading data from file. Function should accept two parameters,
               file name and discharge current. If data contains measured current, discharge current
               could be discarded. Function should return N-by-3 array, where first column contains
               time in seconds, second voltage in volts and third current in amperes.
       [Ur1 Ur2] - Resistance calculation start and end voltage.Voltage is specified as a percentage of
               nominal voltage Ur,e.g.[0.9 0.7]. Straight-line approximation is applied to the discharge
               curve between Ur1 and Ur2. Equivalent serial resistance is calculated from the voltage drop
               at the discharge start time, which is the determined from the value of the straight line
               at the discharge start time.
               Special cases:
                   [1 t] - Ur1 is voltage at first sample after the discharge current was applied;
                           Ur2 is voltage t seconds later.  3rd order polynomial approximation is used to
                           determine the voltage drop at the discharge start time.
       [Uc1 Uc2] - Capacitance calculation start and end voltage.Voltage is specified as a percentage
               of nominal voltage Ur,e.g.[0.9 0.7].
               Special cases:
                   [1 Uc2] - Uc1 is voltage at first sample after the discharge current was applied
                   [1 0] - Uc1 is voltage drop at the discharge start time,Uc2 is voltage of last
                           recorded discharge curve sample
       pf1 - Plot Frequency. Discharge curve plot will be generated each pF1 cycles.
       pf2 - Plot Frequency. Cumulative discharge curve plot will be generated and discharge curve will be
             plotted each pf2 cycles.

    Stored values (in file 'Results.mat'):
       C - Capacitance
       Cs_m - Specific capacitance per mass
       Cs_a - Specific capacitance per area
       R - Equivalent serial resistance
       Rd - Discharge resistance,calculated from self-discharge curve during initial 5s rest period
       nl - Discharge curve nonlinearity,calculated as mean square deviation from the ideal (linear)
            discharge curve.
       Pd - Specific power density, per mass.
       E  - Specific energy, per mass.
    """

    ############# ADD VALIDATION HERE

    noCurrents = len(I) if type(I) is tuple else 1
    noFileGroups = len(files) if type(files) is tuple else 1

    # Determine number of cycles
    f = glob.glob(files[0]) if type(files) is tuple else glob.glob(files)
    noCycles = len(f)
    if noCycles == 0:
        print('No data files selected, check current folder and input number 10, files.')
        return

    # Prealocate result arrays
    C = np.zeros((noCycles,noFileGroups))
    Cs_m = np.zeros((noCycles,noFileGroups))
    Cs_a = np.zeros((noCycles,noFileGroups))
    R = np.zeros((noCycles,noFileGroups))
    Rd = np.zeros((noCycles,noFileGroups))
    nl = np.zeros((noCycles,noFileGroups))
    Pd = np.zeros((noCycles,noFileGroups))
    E = np.zeros((noCycles,noFileGroups))

    # Currents
    for fileIdx in range(noFileGroups):
        if fileIdx < noCurrents:	# if noFileGroups == noCurrents + 1, then measurement is repeated for the first current (repeatability) 
            I0 = I[fileIdx] if type(I) is tuple else I
            fileSuffix = str(I0) + "_mA"
            additionalCurrent = False
        else:   # repeated measurement for first current
            I0 = I[0]
            fileSuffix = str(I0) + "_mA_2nd"
            additionalCurrent = True
        # Replace '.' with '_'
        fileSuffix = fileSuffix.replace('.','_')

        # File list
        f = glob.glob(files[fileIdx]) if type(files) is tuple else glob.glob(files)
        if len(f) != noCycles:
            raise ValueError("All currents should have the same number of cycles (" + str(I0) + " mA: " + str(len(f))
                            + ", should be " + str(noCycles))
        f.sort()    # sorted list of file names

        for n in range(noCycles):
            (c1, r1, rd1, nl1, pd, e) = C_R.C_R(f[n],dataFunct,I0/1000,n+1,noCycles,fileIdx,noCurrents,LimUr,LimUc,pf1,pf2)
            C[n,fileIdx] = c1;
            R[n,fileIdx] = r1;
            Rd[n,fileIdx] = rd1;
            nl[n,fileIdx] = nl1;
            Pd[n,fileIdx] = pd;
            E[n,fileIdx] = e;
        
        if m > 0:
            Cs_m[:,fileIdx] = C[:,fileIdx] / (m / 1000)     # mass in mg
            Pd[:,fileIdx] = Pd[:,fileIdx] / (m / 1000)
            E[:,fileIdx] = E[:,fileIdx] / (m / 1000)
        if S > 0:
            Cs_a[:,fileIdx] = C[:,fileIdx] / S


        baseX = 0.7
        baseY = 0.4
        fontSize = 11

        # CAPACITANCE PLOT
        h = plt.figure(figsize = [10,6], dpi = 96)
        ax = h.add_subplot(1,1,1)
        ax.plot(C[:,fileIdx],'b', linewidth =2)
        bottom, top = plt.ylim()
        plt.ylim(0, 1.05 * top)

        ax.set_title("Capacitance")
        ax.set_xlabel("Cycle")
        ax.set_ylabel("C[F]")
        plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

        annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

        plt.savefig(fname = "_C_" + fileSuffix, dpi = 300)
        if noCurrents > 1:
           plt.close(h)

        # ESR PLOT
        h = plt.figure(figsize = [10,6], dpi = 96)
        ax = h.add_subplot(1,1,1)
        ax.plot(R[:,fileIdx],'b', linewidth =2)
        bottom, top = plt.ylim()
        plt.ylim(0, 1.05 * top)

        ax.set_title("ESR")
        ax.set_xlabel("Cycle")
        ax.set_ylabel("R[\u03A9]")
        plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

        annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

        plt.savefig(fname = "_ESR_" + fileSuffix, dpi = 300)
        if noCurrents > 1:
           plt.close(h)

        # SELF-DISCHARGE RESISTANCE PLOT
        h = plt.figure(figsize = [10,6], dpi = 96)
        ax = h.add_subplot(1,1,1)
        ax.plot(Rd[:,fileIdx],'b', linewidth =2)
        bottom, top = plt.ylim()
        plt.ylim(0, 1.05 * top)

        ax.set_title("Self-discharge resistance")
        ax.set_xlabel("Cycle")
        ax.set_ylabel("Rd[\u03A9]")
        plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

        annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

        plt.savefig(fname = "_Rd_" + fileSuffix, dpi = 300)
        if noCurrents > 1:
           plt.close(h)

        # DISCHARGE CURVE NON-LINEARITY PLOT
        h = plt.figure(figsize = [10,6], dpi = 96)
        ax = h.add_subplot(1,1,1)
        ax.plot(nl[:,fileIdx],'b', linewidth =2)
        bottom, top = plt.ylim()
        plt.ylim(0, 1.05 * top)

        ax.set_title("Discharge curve non-linearity")
        ax.set_xlabel("Cycle")
        ax.set_ylabel("MSE")
        plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

        annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

        plt.savefig(fname = "_Nl_" + fileSuffix, dpi = 300)
        if noCurrents > 1:
           plt.close(h)

        # Cumulative discharge plot by currents
        h = plt.figure(num = C_R.C_R.cumulativeDischargePlot.number)
        plt.legend(loc='upper right')
        plt.savefig(fname = "_discharge_" + fileSuffix, dpi = 300)
        del C_R.C_R.cumulativeDischargePlot
        if noCurrents > 1:
           plt.close(h)
        
        if m > 0:
            # SPECIFIC CAPACITANCE PLOT (by mass)
            h = plt.figure(figsize = [10,6], dpi = 96)
            ax = h.add_subplot(1,1,1)
            ax.plot(Cs_m[:,fileIdx],'b', linewidth =2)
            bottom, top = plt.ylim()
            plt.ylim(0, 1.05 * top)

            ax.set_title("Specific capacitance")
            ax.set_xlabel("Cycle")
            ax.set_ylabel("Cs_m[F/g]")
            plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

            annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

            plt.savefig(fname = "_Cs_m_" + fileSuffix, dpi = 300)
            if noCurrents > 1:
               plt.close(h)

            # MAXIMUM POWER DENSITY PLOT (by mass)
            h = plt.figure(figsize = [10,6], dpi = 96)
            ax = h.add_subplot(1,1,1)
            ax.plot(Pd[:,fileIdx],'b', linewidth =2)
            bottom, top = plt.ylim()
            plt.ylim(0, 1.05 * top)

            ax.set_title("Maximum power density")
            ax.set_xlabel("Cycle")
            ax.set_ylabel("Pdm[W/g]")
            plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

            annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

            plt.savefig(fname = "_Pdm_" + fileSuffix, dpi = 300)
            if noCurrents > 1:
               plt.close(h)

            # ENERGY DENSITY PLOT (by mass)
            h = plt.figure(figsize = [10,6], dpi = 96)
            ax = h.add_subplot(1,1,1)
            ax.plot(E[:,fileIdx],'b', linewidth =2)
            bottom, top = plt.ylim()
            plt.ylim(0, 1.05 * top)

            ax.set_title("Energy density")
            ax.set_xlabel("Cycle")
            ax.set_ylabel("Ed[J/g]")
            plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

            annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

            plt.savefig(fname = "_Ed_" + fileSuffix, dpi = 300)
            if noCurrents > 1:
               plt.close(h)

        if S > 0:
            # SPECIFIC CAPACITANCE PLOT (by area)
            h = plt.figure(figsize = [10,6], dpi = 96)
            ax = h.add_subplot(1,1,1)
            ax.plot(Cs_a[:,fileIdx],'b', linewidth =2)
            bottom, top = plt.ylim()
            plt.ylim(0, 1.05 * top)

            ax.set_title("Specific capacitance")
            ax.set_xlabel("Cycle")
            ax.set_ylabel("Cs_a[F/cm\u00B2]")
            plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

            annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

            plt.savefig(fname = "_Cs_m_" + fileSuffix, dpi = 300)
            if noCurrents > 1:
               plt.close(h)


    # CURRENT INFLUENCE ANALYSIS
    if noCurrents > 1:
        baseX = 0.7
        baseY = 0.4
        fontSize = 11

        # CAPACITANCE PLOT
        h = plt.figure(figsize = [10,6], dpi = 96)
        ax = h.add_subplot(1,1,1)
        ax.plot(I, C[noCycles-1,0:noCurrents],'b', linewidth = 2)
        if additionalCurrent:
            ax.plot(I[0],C[noCycles-1,noFileGroups-1],'rx')
        bottom, top = plt.ylim()
        plt.ylim(0, 1.05 * top)

        ax.set_title("Capacitance")
        ax.set_xlabel("I[mA]")
        ax.set_ylabel("C[F]")
        plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

        annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

        plt.savefig(fname = "_C", dpi = 300)

        # ESR PLOT
        h = plt.figure(figsize = [10,6], dpi = 96)
        ax = h.add_subplot(1,1,1)
        ax.plot(I, R[noCycles-1,0:noCurrents],'b', linewidth =2)
        if additionalCurrent:
            ax.plot(I[0],R[noCycles-1,noFileGroups-1],'rx')
        bottom, top = plt.ylim()
        plt.ylim(0, 1.05 * top)

        ax.set_title("ESR")
        ax.set_xlabel("I[mA]")
        ax.set_ylabel("R[\u03A9]")
        plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

        annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

        plt.savefig(fname = "_ESR", dpi = 300)

        # SELF-DISCHARGE RESISTANCE PLOT
        h = plt.figure(figsize = [10,6], dpi = 96)
        ax = h.add_subplot(1,1,1)
        ax.plot(I, Rd[noCycles-1,0:noCurrents],'b', linewidth =2)
        if additionalCurrent:
            ax.plot(I[0],Rd[noCycles-1,noFileGroups-1],'rx')
        bottom, top = plt.ylim()
        bottom, top = plt.ylim()
        plt.ylim(0, 1.05 * top)

        ax.set_title("Self-discharge resistance")
        ax.set_xlabel("I[mA]")
        ax.set_ylabel("Rd[\u03A9]")
        plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

        annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

        plt.savefig(fname = "_Rd", dpi = 300)

        # DISCHARGE CURVE NON-LINEARITY PLOT
        h = plt.figure(figsize = [10,6], dpi = 96)
        ax = h.add_subplot(1,1,1)
        ax.plot(I, nl[noCycles-1,0:noCurrents],'b', linewidth =2)
        if additionalCurrent:
            ax.plot(I[0],nl[noCycles-1,noFileGroups-1],'rx')
        bottom, top = plt.ylim()
        plt.ylim(0, 1.05 * top)

        ax.set_title("Discharge curve non-linearity")
        ax.set_xlabel("I[mA]")
        ax.set_ylabel("MSE")
        plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

        annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

        plt.savefig(fname = "_Nl", dpi = 300)

        
        if m > 0:
            # SPECIFIC CAPACITANCE PLOT (by mass)
            h = plt.figure(figsize = [10,6], dpi = 96)
            ax = h.add_subplot(1,1,1)
            ax.plot(I, Cs_m[noCycles-1,0:noCurrents],'b', linewidth =2)
            if additionalCurrent:
                ax.plot(I[0],Cs_m[noCycles-1,noFileGroups-1],'rx')
            bottom, top = plt.ylim()
            plt.ylim(0, 1.05 * top)

            ax.set_title("Specific capacitance")
            ax.set_xlabel("I[mA]")
            ax.set_ylabel("Cs_m[F/g]")
            plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

            annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

            plt.savefig(fname = "_Cs_m", dpi = 300)

            # MAXIMUM POWER DENSITY PLOT (by mass)
            h = plt.figure(figsize = [10,6], dpi = 96)
            ax = h.add_subplot(1,1,1)
            ax.plot(I, Pd[noCycles-1,0:noCurrents],'b', linewidth =2)
            if additionalCurrent:
                ax.plot(I[0],Pd[noCycles-1,noFileGroups-1],'rx')
            bottom, top = plt.ylim()
            plt.ylim(0, 1.05 * top)

            ax.set_title("Maximum power density")
            ax.set_xlabel("I[mA]")
            ax.set_ylabel("Pdm[W/g]")
            plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

            annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

            plt.savefig(fname = "_Pdm", dpi = 300)
            if noCurrents > 1:
               plt.close(h)

            # ENERGY DENSITY PLOT (by mass)
            h = plt.figure(figsize = [10,6], dpi = 96)
            ax = h.add_subplot(1,1,1)
            ax.plot(I,E[noCycles-1,0:noCurrents],'b', linewidth =2)
            if additionalCurrent:
                ax.plot(I[0],E[noCycles-1,noFileGroups-1],'rx')
            bottom, top = plt.ylim()
            plt.ylim(0, 1.05 * top)

            ax.set_title("Energy density")
            ax.set_xlabel("I[mA]")
            ax.set_ylabel("Ed[J/g]")
            plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

            annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

            plt.savefig(fname = "_Ed", dpi = 300)
            if noCurrents > 1:
               plt.close(h)

        if S > 0:
            # SPECIFIC CAPACITANCE PLOT (by area)
            h = plt.figure(figsize = [10,6], dpi = 96)
            ax = h.add_subplot(1,1,1)
            ax.plot(I, Cs_a[noCycles-1,0:noCurrents],'b', linewidth =2)
            if additionalCurrent:
                ax.plot(I[0],Cs_a[noCycles-1,noFileGroups-1],'rx')
            bottom, top = plt.ylim()
            plt.ylim(0, 1.05 * top)

            ax.set_title("Specific capacitance")
            ax.set_xlabel("I[mA]")
            ax.set_ylabel("Cs_a[F/cm\u00B2]")
            plt.grid(axis="both", color=(0.25, 0.25, 0.25), linestyle=':')

            annotate.annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);

            plt.savefig(fname = "_Cs_m", dpi = 300)
            if noCurrents > 1:
               plt.close(h)

        # Cumulative discharge plot
        plt.figure(num = C_R.C_R.dischargeIplot.number)
        plt.legend(loc='upper right')
        plt.savefig(fname = "_discharge_I_", dpi = 300)


    plt.show()      # display all figures
        

# SAMPLE CALLS
#os.chdir(r"D:\Projekti\FKIT\Superkondenzator\Testiranje\Suzana mjerenja\1012")
#SCanalysis(2.7,22.1,6.9,2,"label1","label2","label3","label4","label5", "1012KRONO2_*.oxw", \
#    "readFunctions.readAL",(1,1),(0.9,0.7),10,20)
#os.chdir(r"D:\Projekti\FKIT\Superkondenzator\Testiranje\Suzana mjerenja\2017-02 Ljubljana")
#SCanalysis(2.7,20,27.7,2,"80% AC, 10% CB, 10% PVDF", "80oC / 5MPa press", "115um thickness", \
#    "packed in glove box", "", "*KRONO2*.oxw", "readFunctions.readAL", (0.9,0.7), (0.9,0.7), 10, 50)
#os.chdir(r"D:\Projekti\FKIT\Superkondenzator\Testiranje\Suzana mjerenja\Strujni testovi")
#SCanalysis(2.7, (8, 9.8, 12, 14.7, 18.1, 22.1, 27.1, 33.3, 40.8, 50), 10, 2, \
#    "Line 1", "Line 2", "Line 3", "Line 4", "Line 5", \
#    ("008_0mA_*.oxw", "009_8mA_*.oxw", "012_0mA_*.oxw", "014_7mA_*.oxw", "018_1mA_*.oxw", \
#    "022_1mA_*.oxw", "027_1mA_*.oxw", "033_3mA_*.oxw", "040_8mA_*.oxw", "050_0mA_*.oxw", "nakon_*.oxw"), \
#    "readFunctions.readAL", (1,1), (0.9,0.7), 1, 1)
