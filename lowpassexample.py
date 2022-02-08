import numpy as np
import matplotlib.pyplot as plt
np.random.seed(123456789)

tspan =  np.linspace(0, 10*np.pi, 1000) 

Vin = lambda t: np.sin(t)
mynoise = 1*np.random.randn(*Vin(tspan).shape)
Vinarray = Vin(tspan) + mynoise

def spectrum(x):
    return np.abs(np.fft.ifft(x))

plt.plot(tspan, Vinarray)
plt.plot(tspan, Vin(tspan))

mysignal_spectrum = spectrum(Vin(tspan))
mynoise_spectrum  = spectrum(mynoise)
N1 = 500
f = np.arange(N1)
plt.figure(figsize=(8,4))
plt.plot(mynoise_spectrum[:N1])
plt.plot(f,mysignal_spectrum[:N1])
plt.xlim(0,N1)
plt.legend(('signal','noise'))
plt.xlabel('frequency')
plt.ylabel('amplitude')
plt.show()


Voutarray = np.zeros_like(Vinarray)

def lowpass(tstep, Vin, tau):
    Vout = np.zeros_like(Vin)
    for t in range(Vin.size-1):
        Vout[t+1] = Vout[t]+(tstep/tau)*(Vin[t]-Vout[t])
    return Vout
    

tau = 1
tstep = tspan[1]-tspan[0]
plt.plot(tspan,lowpass(tstep, Vinarray, tau))
plt.show()




