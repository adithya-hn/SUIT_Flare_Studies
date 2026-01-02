import pywt
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import morlet2
print(pywt.wavelist(family=None, kind='continuous'))

w = pywt.Wavelet('db5')
phi, psi, x = w.wavefun(level=6)

plt.plot(x, psi,'ko-')
plt.title("Daubechies db4 wavelet (ψ)")
plt.xlabel("t")
plt.ylabel("Amplitude")
plt.show()


# import numpy as np
# from scipy.signal import morlet2

# t = np.linspace(-5, 5, 1000)
# psi = morlet2(M=len(t), s=1.0, w=6)

# plt.plot(t, np.real(psi), label="Real")
# plt.plot(t, np.imag(psi), '--', label="Imag")
# plt.legend()
# plt.title("Morlet wavelet")
# plt.show()
"""

t = np.arange(200)
signal = np.zeros_like(t, dtype=float)
signal[80:85] = 1.0  # impulsive brightening

plt.plot(signal)
plt.title("Toy transient signal")
plt.show()


coeffs = pywt.wavedec(signal, 'db4', level=4)

for i in range(1, 5):
    plt.plot(coeffs[-i], label=f"D{i}")
plt.legend()
plt.title("DWT detail coefficients")
plt.show()


from scipy.signal import cwt

scales = np.arange(1, 64)
cwt_coeff = cwt(signal, morlet2, scales, w=6)

plt.imshow(np.abs(cwt_coeff)**2,
           aspect='auto',
           cmap='inferno',
           origin='lower')
plt.ylabel("Scale")
plt.xlabel("Time")
plt.title("CWT scalogram")
plt.colorbar(label="Power")
plt.show()
"""