import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

YR   = 365.2563

name = "3_bodies.txt"

try:
    data = np.genfromtxt(name)
except ValueError as err:
    try:
        data = np.genfromtxt(name, skip_footer=1)
        print("Still writing...")
    except:
        raise ValueError (err)
print("Shape:", data.shape)
N = data.shape[0]

df = pd.DataFrame(data)
df = df.rename(columns=dict(
    zip(df.columns,
     ["s0", "o0",
     "a1", "K1", "s1", "o1", "H1",
     "a2", "K2", "s2", "o2", "H2",
     "n1", "n2",
     "t", "dt", "dt_adap"])
    ))
df.t       /= YR
df.dt      /= YR
df.dt_adap /= YR
df["e1"] = np.sqrt(df.K1**2 + df.H1**2)
df["e2"] = np.sqrt(df.K2**2 + df.H2**2)
df["w1"] = np.arctan2(df.H1, df.K1)
df["w2"] = np.arctan2(df.H2, df.K2)

name2 = "2_bodies.txt"
try:
    data2 = np.genfromtxt(name2)
except ValueError as err:
    try:
        data2 = np.genfromtxt(name2, skip_footer=1)
        print("Still writing...")
    except:
        raise ValueError (err)
print("Shape:", data2.shape)
N2 = data2.shape[0]

df2 = pd.DataFrame(data2)
df2 = df2.rename(columns=dict(
    zip(df2.columns,
     ["s0", "o0",
     "a1", "K1", "s1", "o1", "H1",
     "a2", "K2", "s2", "o2", "H2",
     "n1", "n2",
     "t", "dt", "dt_adap"])
    ))
df2.t       /= YR
df2.dt      /= YR
df2.dt_adap /= YR
df2["e1"] = np.sqrt(df2.K1**2 + df2.H1**2)
df2["e2"] = np.sqrt(df2.K2**2 + df2.H2**2)
df2["w1"] = np.arctan2(df2.H1, df2.K1)
df2["w2"] = np.arctan2(df2.H2, df2.K2)

plt.figure(dpi=150)
plt.title("TimeStep")
plt.plot(df.t, df.dt, '.', ms=0.5, label='dt')
plt.plot(df.t, df.dt_adap, '.', ms=0.5, label='dt_adap')
plt.xlabel("Time / yr")
plt.ylabel("dt / yr")
plt.semilogx()
plt.semilogy()
plt.legend()
plt.tight_layout()

plt.figure(dpi=200, figsize=(8,3))
plt.subplot(1,2,1)
plt.title("a$_1$")
plt.plot(df.t, df.a1, '.', ms=0.5,  label="3 bodies")
plt.plot(df2.t, df2.a1, '.', ms=0.5,  label="2 bodies")
plt.xlabel("Time / yr")
plt.ylabel("a / AU")
plt.semilogx()
plt.legend()
plt.subplot(1,2,2)
plt.title("a$_2$")
plt.plot(df.t, df.a2, '.', ms=0.5,  label='3 bodies')
plt.plot(df2.t, df2.a2, '.', ms=0.5,  label='2 bodies')
plt.xlabel("Time / yr")
plt.semilogx()
plt.legend()
plt.tight_layout()

plt.figure(dpi=150)
plt.title("a$_1$")
plt.plot(df.t * 1e-9, df.a1, '.', ms=0.5,  label="3 bodies")
plt.plot(df2.t * 1e-9, df2.a1, '.', ms=0.5,  label="2 bodies")
plt.xlabel("Time / Gyr")
plt.ylabel("a / AU")
plt.xlim(0, 1)
plt.ylim(0.0475, 0.05)
plt.legend()
plt.tight_layout()

plt.figure(dpi=150)
plt.title("a$_1$")
plt.plot(df.t * 1e-9, df.a1, '.', ms=0.5,  label="3 bodies")
plt.plot(df2.t * 1e-9, df2.a1, '.', ms=0.5,  label="2 bodies")
plt.xlabel("Time / Gyr")
plt.ylabel("a / AU")
plt.xlim(0, 5)
plt.legend()
plt.tight_layout()

plt.figure(dpi=200, figsize=(8,3))
plt.subplot(1,2,1)
plt.title("e$_1$")
plt.plot(df.t, df.e1, '.', ms=0.5,  label="3 bodies")
plt.plot(df2.t, df2.e1, '.', ms=0.5,  label="2 bodies")
plt.xlabel("Time / yr")
plt.ylabel("e")
plt.semilogx()
plt.legend()
plt.subplot(1,2,2)
plt.title("e$_2$")
plt.plot(df.t, df.e2, '.', ms=0.5,  label="3 bodies")
plt.plot(df2.t, df2.e2, '.', ms=0.5,  label="2 bodies")
plt.xlabel("Time / yr")
plt.semilogx()
plt.legend()
plt.tight_layout()

plt.figure(dpi=150)
plt.plot(df.t * 1e-9, df.e1, '.', ms=0.5, label="e$_1$ - 3 bodies")
plt.plot(df2.t * 1e-9, df2.e1, '.', ms=0.5, label="e$_1$ - 2 bodies")
plt.plot(df.t * 1e-9, df.e2, '.', ms=0.5, label="e$_2$ - 3 bodies")
plt.plot(df2.t * 1e-9, df2.e2, '.', ms=0.5, label="e$_2$ - 2 bodies")
plt.xlabel("Time / Gyr")
plt.ylabel("e")
plt.xlim(0, 0.5)
plt.ylim(0, 0.12)
plt.legend()
plt.tight_layout()

plt.figure(dpi=200, figsize=(10,3))
plt.subplot(1,3,1)
plt.title("$\Omega_0 / n_1$")
plt.plot(df.t, df.s0/df.n1, '.', ms=0.5, label="3 bodies")
plt.plot(df2.t, df2.s0/df2.n1, '.', ms=0.5, label="2 bodies")
plt.xlabel("Time / yr")
plt.ylabel("$\Omega / n$")
plt.semilogx()
plt.legend()
plt.subplot(1,3,2)
plt.title("$\Omega_1 / n_1$")
plt.plot(df.t, df.s1/df.n1, '.', ms=0.5, label="3 bodies")
plt.plot(df2.t, df2.s1/df2.n1, '.', ms=0.5, label="2 bodies")
plt.xlabel("Time / yr")
plt.semilogx()
plt.semilogy()
plt.legend()
plt.subplot(1,3,3)
plt.title("$\Omega_2 / n_2$")
plt.plot(df.t, df.s2/df.n2, '.', ms=0.5, label="3 bodies")
plt.plot(df2.t, df2.s2/df2.n2, '.', ms=0.5, label="2 bodies")
plt.xlabel("Time / yr")
plt.semilogx()
plt.semilogy()
plt.legend()
plt.tight_layout()

plt.figure(dpi=200, figsize=(10,3))
plt.subplot(1,3,1)
plt.title("$\epsilon_0$")
plt.plot(df.t, df.o0, '.', ms=0.5, label="3 bodies")
plt.plot(df2.t, df2.o0, '.', ms=0.5, label="2 bodies")
plt.xlabel("Time / yr")
plt.ylabel("$\epsilon$")
plt.semilogx()
plt.subplot(1,3,2)
plt.title("$\epsilon_1$")
plt.plot(df.t, df.o1, '.', ms=0.5, label="3 bodies")
plt.plot(df2.t, df2.o1, '.', ms=0.5, label="2 bodies")
plt.xlabel("Time / yr")
plt.semilogx()
plt.semilogy()
plt.legend()
plt.subplot(1,3,3)
plt.title("$\epsilon_2$")
plt.plot(df.t, df.o2, '.', ms=0.5, label="3 bodies")
plt.plot(df2.t, df2.o2, '.', ms=0.5, label="2 bodies")
plt.xlabel("Time / yr")
plt.semilogx()
plt.semilogy()
plt.legend()
plt.tight_layout()

plt.figure(dpi=200, figsize=(8,3))
plt.subplot(1,2,1)
plt.title("$\omega_1$")
plt.plot(df.t, df.w1, '.', ms=0.5, label="3 bodies")
plt.plot(df2.t, df2.w1, '.', ms=0.5, label="2 bodies")
plt.xlabel("Time / yr")
plt.ylabel("$\omega$")
plt.semilogx()
plt.legend()
plt.subplot(1,2,2)
plt.title("$\omega_2$")
plt.plot(df.t, df.w2, '.', ms=0.5, label="3 bodies")
plt.plot(df2.t, df2.w2, '.', ms=0.5, label="2 bodies")
plt.xlabel("Time / yr")
plt.semilogx()
plt.legend()
plt.tight_layout()

plt.figure(dpi=200, figsize=(8,3))
plt.subplot(1,2,1)
plt.title("$H_1,K_1$")
plt.plot(df.K1, df.H1, '.', ms=0.5, label="3 bodies")
plt.plot(df2.K1, df2.H1, '.', ms=0.5, label="2 bodies")
plt.xlabel("K1")
plt.ylabel("H1")
plt.axis("equal")
plt.legend()
plt.subplot(1,2,2)
plt.title("$H_2,K_2$")
plt.plot(df.K2, df.H2, '.', ms=0.5, label="3 bodies")
plt.plot(df2.K2, df2.H2, '.', ms=0.5, label="2 bodies")
plt.xlabel("K2")
plt.ylabel("H2")
plt.axis("equal")
plt.legend()
plt.tight_layout()

plt.show()
