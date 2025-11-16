#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    ex1rk4 = pd.read_csv("rk4.csv", sep=",", header=None);
    t_rk4 = ex1rk4.iloc[:,0].to_numpy();
    y_rk4 = ex1rk4.iloc[:,1].to_numpy();
    ydot_rk4 = ex1rk4.iloc[:,2].to_numpy();
    ex1euler = pd.read_csv("euler.csv", sep=",", header=None);
    t_euler = ex1euler.iloc[:,0].to_numpy();
    y_euler = ex1euler.iloc[:,1].to_numpy();
    # ydot_euler = ex1euler.iloc[:,2].to_numpy();
    t = np.linspace(0,6,10000);
    y = t**3 + t**2;
    fig1 = plt.figure(layout="constrained");
    plt.grid(alpha=0.4);
    plt.title("Comparação entre RK4, Euler, e solução analítica")
    plt.xlim([0,6]);
    plt.ylim([0,252]);
    plt.xlabel(r"$t$");
    plt.ylabel(r"$y$");
    plt.plot(t, y, label=r"$y(t) = t^3 + t^2$", color="#ca9ee6");
    plt.plot(t_rk4, y_rk4, label="RK4", color="#ea76cb");
    plt.plot(t_euler, y_euler, label="Euler", color="#8839ef");
    plt.legend();
    plt.show();
    fig1.savefig("1.png", dpi=500);
    fig1b = plt.figure(layout="constrained");
    plt.title("Erro relativo entre solução analítica e RK4")
    plt.xlim([0,6]);
    plt.xlabel(r"$t$");
    plt.ylabel(r"$\varepsilon$");
    plt.plot(t_rk4, 1 - y_rk4/(t_rk4**3 + t_rk4**2), label=r"$\varepsilon_y(t) = \frac{y(t) - y_{RK4}(t)}{y(t)}$", color="#ea76cb");
    plt.plot(t_rk4, 1 - ydot_rk4/(3*t_rk4**2 + 2*t_rk4), label=r"$\varepsilon_{\dot{y}}(t) = \frac{\dot{y}(t) - \dot{y}_{RK4}(t)}{\dot{y}(t)}$", color="#8839ef");
    plt.legend();
    plt.show();
    fig1b.savefig("1b.png", dpi=500);

    ex2_a01 = pd.read_csv("2a01.csv", sep=",", header = None);
    x_a01 = ex2_a01.iloc[:,1].to_numpy();
    v_a01 = ex2_a01.iloc[:,2].to_numpy();
    ex2_a05 = pd.read_csv("2a05.csv", sep=",", header = None);
    x_a05 = ex2_a05.iloc[:,1].to_numpy();
    v_a05 = ex2_a05.iloc[:,2].to_numpy();
    ex2_a10 = pd.read_csv("2a10.csv", sep=",", header = None);
    x_a10 = ex2_a10.iloc[:,1].to_numpy();
    v_a10 = ex2_a10.iloc[:,2].to_numpy();
    fig2a = plt.figure(layout="constrained");
    plt.grid(alpha=0.4);
    plt.title("Espaço de fase para o potencial poço duplo")
    plt.xlabel(r"$x$");
    plt.ylabel(r"$v$");
    plt.plot(x_a01, v_a01, label=r"$\dot{x}(0) = 0.1$", color="#ca9ee6");
    plt.plot(x_a05, v_a05, label=r"$\dot{x}(0) = 0.5$", color="#ea76cb");
    plt.plot(x_a10, v_a10, label=r"$\dot{x}(0) = 1.0$", color="#8839ef");
    plt.legend();
    plt.show();
    fig2a.savefig("2a.png", dpi=500);

    ex2_b25 = pd.read_csv("2b25.csv", sep=",", header = None);
    x_b25 = ex2_b25.iloc[:,1].to_numpy();
    v_b25 = ex2_b25.iloc[:,2].to_numpy();
    ex2_b80 = pd.read_csv("2b80.csv", sep=",", header = None);
    x_b80 = ex2_b80.iloc[:,1].to_numpy();
    v_b80 = ex2_b80.iloc[:,2].to_numpy();
    fig2b = plt.figure(layout="constrained");
    plt.grid(alpha=0.4);
    plt.title("Espaço de fase para o potencial poço duplo amortecido")
    plt.xlabel(r"$x$");
    plt.ylabel(r"$v$");
    plt.plot(x_b25, v_b25, label=r"$2\gamma = 0.25$", color="#ca9ee6");
    plt.plot(x_b80, v_b80, label=r"$2\gamma = 0.80$", color="#ea76cb");
    plt.legend();
    plt.show();
    fig2b.savefig("2b.png", dpi=500);

    ex2_c190 = pd.read_csv("2c190.csv", sep=",", header = None);
    x_c190 = ex2_c190.iloc[:,1].to_numpy();
    v_c190 = ex2_c190.iloc[:,2].to_numpy();
    ex2_c203 = pd.read_csv("2c203.csv", sep=",", header = None);
    x_c203 = ex2_c203.iloc[:,1].to_numpy();
    v_c203 = ex2_c203.iloc[:,2].to_numpy();
    ex2_c240 = pd.read_csv("2c240.csv", sep=",", header = None);
    x_c240 = ex2_c240.iloc[:,1].to_numpy();
    v_c240 = ex2_c240.iloc[:,2].to_numpy();
    ex2_c330 = pd.read_csv("2c330.csv", sep=",", header = None);
    x_c330 = ex2_c330.iloc[:,1].to_numpy();
    v_c330 = ex2_c330.iloc[:,2].to_numpy();
    ex2_c600 = pd.read_csv("2c600.csv", sep=",", header = None);
    x_c600 = ex2_c600.iloc[:,1].to_numpy();
    v_c600 = ex2_c600.iloc[:,2].to_numpy();
    fig2c = plt.figure(layout="constrained");
    plt.grid(alpha=0.4);
    plt.title("Espaço de fase para o potencial poço duplo amortecido forçado")
    plt.xlabel(r"$x$");
    plt.ylabel(r"$v$");
    plt.plot(x_c190, v_c190, label=r"$F = 0.190$", color = "#d20f39");
    plt.plot(x_c203, v_c203, label=r"$F = 0.203$", color = "#8839ef");
    plt.plot(x_c240, v_c240, label=r"$F = 0.240$", color = "#40a02b");
    plt.plot(x_c330, v_c330, label=r"$F = 0.330$", color = "#fe640b");
    plt.plot(x_c600, v_c600, label=r"$F = 0.600$", color = "#7287fd");
    plt.legend();
    plt.show();
    fig2c.savefig("2c.png", dpi=500);

    poincaré_section = pd.read_csv("poincaré_section.csv", sep=",", header=None);
    f = poincaré_section.iloc[:,3].to_numpy();
    x = poincaré_section.iloc[:,1].to_numpy();
    fig2d = plt.figure(layout="constrained");
    plt.title("Diagrama de bifurcação")
    plt.xlabel(r"$F$");
    plt.ylabel(r"$x$");
    plt.scatter(f, x, s = 0.1, color="#ea76cb");
    plt.show();
    fig2d.savefig("2d.png", dpi=500);

    poincaré_map = pd.read_csv("poincaré_map.csv", sep=",", header=None);
    x = poincaré_map.iloc[:,1].to_numpy();
    v = poincaré_map.iloc[:,2].to_numpy();
    fig2e = plt.figure(layout="constrained");
    plt.title("Mapa de Poincaré")
    plt.xlabel(r"$x$");
    plt.ylabel(r"$v$");
    plt.scatter(x, v, s = 0.1, color="#ea76cb");
    plt.show();
    fig2e.savefig("2e.png", dpi=500);
