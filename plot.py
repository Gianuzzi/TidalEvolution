import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

YR = 365.2563

#name = input("Filename: ")
name = "Salida.txt"

try:
    data = np.genfromtxt(name)
except ValueError as err:
    try:
        data = np.genfromtxt(name, skip_footer=1)
        print("Still writing...")
    except:
        raise ValueError (err)
print("Shape:", data.shape)

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
df["e1"]    = np.sqrt(df.K1**2 + df.H1**2)
df["e2"]    = np.sqrt(df.K2**2 + df.H2**2)
df["w1"]    = np.arctan2(df.H1, df.K1)
df["w2"]    = np.arctan2(df.H2, df.K2)

plt.figure(dpi=100)
plt.plot(df.t, df.dt, ',', label='dt')
plt.plot(df.t, df.dt_adap, ',', label='dt_adap')
plt.xlabel("Time/yr")
plt.ylabel("dt/yr")
plt.semilogx()
plt.semilogy()
plt.legend()
plt.tight_layout()
# plt.show()


plt.figure(dpi=150, figsize=(8,3))
plt.subplot(1,2,1)
plt.title("$m_1$")
plt.plot(df.t, df.a1, ',')
plt.xlabel("Time/yr")
plt.ylabel("a / AU")
plt.semilogx()
plt.subplot(1,2,2)
plt.title("$m_2$")
plt.plot(df.t, df.a2, ',')
plt.xlabel("Time/yr")
plt.semilogx()
plt.tight_layout()
# plt.show()

plt.figure(dpi=150, figsize=(8,3))
plt.subplot(1,2,1)
plt.title("$m_1$")
plt.plot(df.t, df.e1, ',')
plt.xlabel("Time/yr")
plt.ylabel("e")
plt.semilogx()
plt.subplot(1,2,2)
plt.title("$m_2$")
plt.plot(df.t, df.e2, ',')
plt.xlabel("Time/yr")
plt.semilogx()
plt.tight_layout()
# plt.show()

plt.figure(dpi=150, figsize=(10,3))
plt.subplot(1,3,1)
plt.title("$m_0$")
plt.plot(df.t, df.s0/df.n1, ',')
plt.xlabel("Time/yr")
plt.ylabel("$\Omega / n$")
plt.semilogx()
plt.subplot(1,3,2)
plt.title("m1")
plt.plot(df.t, df.s1/df.n1, ',')
plt.xlabel("Time/yr")
plt.semilogx()
plt.semilogy()
plt.subplot(1,3,3)
plt.title("m2")
plt.plot(df.t, df.s2/df.n2, ',')
plt.xlabel("Time/yr")
plt.semilogx()
plt.semilogy()
plt.tight_layout()
# plt.show()

plt.figure(dpi=150, figsize=(10,3))
plt.subplot(1,3,1)
plt.title("$m_0$")
plt.plot(df.t, df.o0, ',')
plt.xlabel("Time/yr")
plt.ylabel("$\epsilon$")
plt.semilogx()
plt.subplot(1,3,2)
plt.title("$m_1$")
plt.plot(df.t, df.o1, ',')
plt.xlabel("Time/yr")
plt.semilogx()
plt.semilogy()
plt.subplot(1,3,3)
plt.title("$m_2$")
plt.plot(df.t, df.o2, ',')
plt.xlabel("Time/yr")
plt.semilogx()
plt.semilogy()
plt.tight_layout()
# plt.show()

plt.figure(dpi=150, figsize=(8,3))
plt.subplot(1,2,1)
plt.title("$m_1$")
plt.plot(df.t, df.w1, ',')
plt.xlabel("Time/yr")
plt.ylabel("$\omega$")
plt.semilogx()
plt.subplot(1,2,2)
plt.title("$m_2$")
plt.plot(df.t, df.w2, ',')
plt.xlabel("Time/yr")
plt.semilogx()
plt.tight_layout()
# plt.show()

plt.figure(dpi=150, figsize=(8,3))
plt.subplot(1,2,1)
plt.title("$m_1$")
plt.plot(df.t, df.K1, ',')
plt.semilogx()
plt.xlabel("Time/yr")
plt.ylabel("K")
plt.subplot(1,2,2)
plt.title("$m_2$")
plt.plot(df.t, df.K2, ',')
plt.semilogx()
plt.xlabel("Time/yr")
plt.tight_layout()
# plt.show()

plt.figure(dpi=150, figsize=(8,3))
plt.subplot(1,2,1)
plt.title("$m_1$")
plt.plot(df.t, df.H1, ',')
plt.semilogx()
plt.xlabel("Time/yr")
plt.ylabel("H")
plt.subplot(1,2,2)
plt.title("$m_2$")
plt.plot(df.t, df.H2, ',')
plt.semilogx()
plt.xlabel("Time/yr")
plt.tight_layout()
# plt.show()

plt.figure(dpi=150, figsize=(8,3))
plt.subplot(1,2,1)
plt.title("$m_1$")
plt.plot(df.K1, df.H1, ',')
plt.xlabel("K1")
plt.ylabel("H1")
plt.axis("equal")
plt.subplot(1,2,2)
plt.title("$m_2$")
plt.plot(df.K2, df.H2, ',')
plt.xlabel("K2")
plt.ylabel("H2")
plt.axis("equal")
plt.tight_layout()
plt.show()