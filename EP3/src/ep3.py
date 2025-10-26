#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    ex1 = pd.read_csv("1.csv", sep=",", header=None);
    p = ex1.iloc[:,0].to_numpy();
    e_float = ex1.iloc[:,4].to_numpy();
    e_double = ex1.iloc[:,7].to_numpy();
    fig1 = plt.figure(layout="constrained");
    f_simpson_lin = np.polyfit(p[0:4],e_float[0:4],1);
    d_simpson_lin = np.polyfit(p[0:10],e_double[0:10],1);
    p_f_simpson = np.linspace(0,5);
    p_d_simpson = np.linspace(0,12);
    f_simpson_fit = f_simpson_lin[0]*p_f_simpson + f_simpson_lin[1];
    d_simpson_fit = d_simpson_lin[0]*p_d_simpson + d_simpson_lin[1];
    f_roundoff_lin = np.polyfit(p[4:19],e_float[4:19],1);
    d_roundoff_lin = np.polyfit(p[11:],e_double[11:],1);
    p_float_roundoff = np.linspace(4,21);
    p_double_roundoff = np.linspace(11,26);
    f_roundoff_fit = f_roundoff_lin[0]*p_float_roundoff + f_roundoff_lin[1];
    d_roundoff_fit = d_roundoff_lin[0]*p_double_roundoff + d_roundoff_lin[1];
    plt.grid(alpha=0.4);
    plt.title("Erro no método de Simpson")
    plt.xlim([0,26]);
    plt.ylim([-51,1]);
    plt.xlabel("p");
    plt.ylabel(r"$\log_2(\epsilon)$");
    plt.plot(p_d_simpson, d_simpson_fit, color="#ca9ee6aa", lw=0.5, 
             label=rf"$\Delta$y={d_simpson_lin[0]:.3f}$\Delta$x");
    plt.plot(p_double_roundoff, d_roundoff_fit, color="#ca9ee6aa", ls="-.",
             lw=0.5,label=rf"$\Delta$y={d_roundoff_lin[0]:.3f}$\Delta$x");
    plt.plot(p_float_roundoff, f_roundoff_fit, color="#f4b8e4aa", lw=0.5, 
             label=rf"$\Delta$y={f_roundoff_lin[0]:.3f}$\Delta$x");
    plt.scatter(p, e_float, label="float", color="#ea76cb", marker=".");
    plt.scatter(p, e_double, label="double", color="#8839ef", marker=".");
    plt.legend();
    fig1.savefig("1.png", dpi=500);

    ex2 = pd.read_csv("2.csv", sep=",", header = None);
    t = ex2.iloc[:,1].to_numpy();
    r = ex2.iloc[:,2].to_numpy();
    fig2 = plt.figure(layout="constrained");
    plt.grid(alpha=0.4);
    plt.title("Período de um pêndulo simples")
    plt.xlim([0,180]);
    plt.xlabel(r"$\theta_0$ (°)");
    plt.ylabel(r"$T/T_{\mathrm{Galileu}}$");
    plt.plot(t, r, color="#ca9ee6aa", lw=0.5);
    plt.scatter(t, r, color="#8839ef", marker=".");
    fig2.savefig("2.png",dpi=500);

    ex3 = pd.read_csv("3.csv", sep=",", header=None);
    n = ex3.iloc[:,0].to_numpy();
    i = ex3.iloc[:,1].to_numpy();
    s = ex3.iloc[:,3].to_numpy();
    fig3 = plt.figure(layout="constrained");
    plt.grid(alpha=0.4);
    plt.title(r"Quadratura da parábola $y = 1 - x^2$")
    plt.xlabel(r"$N_t$");
    plt.ylabel(r"$I_m$");
    plt.yticks([18/30, 19/30, 20/30, 21/30],
               labels=["0.6", "", r"$\frac{2}{3}$", "0.7"]);
    plt.ylim([18/30, 21/30]);
    plt.errorbar(x=n, y=i, yerr=s, color="#ca9ee6aa", marker=".",
                 markerfacecolor="#8839ef", linewidth=0.5, elinewidth=1,
                 ecolor="#6c6f85");
    plt.xscale("log",base=2);
    fig3.savefig("3.png",dpi=500);
