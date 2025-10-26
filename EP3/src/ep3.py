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
    plt.grid(alpha=0.4);
    plt.title("Erro no método de Simpson")
    plt.xlim([0,26]);
    plt.xlabel("p");
    plt.ylabel(r"$\log_2(\epsilon)$");
    plt.plot(p, e_float, color="#f4b8e4aa", lw=0.5);
    plt.plot(p, e_double, color="#ca9ee6aa", lw=0.5);
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
    plt.errorbar(x=n, y=i, yerr=s, color="#ca9ee6aa", marker=".", markerfacecolor="#8839ef", linewidth=0.5, elinewidth=1, ecolor="#6c6f85");
    plt.xscale("log",base=2);
    fig3.savefig("3.png",dpi=500);
