import numpy as np
import numpy.random as nram 
import matplotlib.pyplot as plt

# Number of given subcarriers 
nSubcarriers=64
nblocks=10000

# Variance
No = 1

Eb_db=np.arange(1,61,6)
Eb= 10**(Eb_db/10)
SNR=2*(Eb/No)
SNR_db=10*np.log10(SNR)
BER=np.zeros(len(Eb_db))
BER_theory=np.zeros(len(Eb_db))

# Number of unit gain channel taps
L=3

# Number of symbols used for cyclic prefix 
L_tilde=6

for i in range(nblocks):
    # Generation of inphase bits and quadrature bits
    inphase_bits = nram.randint(2,size=nSubcarriers)
    quad_bits = nram.randint(2,size=nSubcarriers)

    # Generation of QPSK symbols
    Qpsk_sym = (2*inphase_bits-1)+1j*(2*quad_bits-1)

    # Generation of noise for both inphase and quadrature bits
    inphase_noise = nram.normal(0,np.sqrt(No/2),nSubcarriers+L-1+L_tilde)
    quad_noise = nram.normal(0,np.sqrt(No/2),nSubcarriers+L-1+L_tilde)

    

    # Complex noise
    noise = inphase_noise+1j*quad_noise

    # Unit gain channel taps
    h=nram.normal(0,np.sqrt(1/2),L)+1j*nram.normal(0,np.sqrt(1/2),L)

    # FFT of unit gain channel taps
    fft_of_h = np.fft.fft(h,nSubcarriers)

    for j in range(len(Eb_db)):
        # Tramitted symbols
        x_of_n=np.sqrt(Eb[j])*Qpsk_sym

        # IFFT of transmitted symbols
        tx_sym=np.fft.ifft(x_of_n)

# Copying the last L_tilde symbols from tansmitted symbols for cyclic   prefix
        cyclic_pref=tx_sym[nSubcarriers-L_tilde:]

        # Addition of cyclic prefix
        tx_sym_cyclic_pref=np.concatenate((cyclic_pref,tx_sym))

        # Received symbols
        rx_sym_cyclic_pref=np.convolve(tx_sym_cyclic_pref,h)+noise

        # Removal of cyclic prefix
        rx_sym=rx_sym_cyclic_pref[L_tilde:L_tilde+nSubcarriers]

        # FFT of received symbols 
        fft_of_rx_sym=np.fft.fft(rx_sym)

        # Single tap equalisation
        sin_tap_equalisation=1/fft_of_h*fft_of_rx_sym

        # Detection of inphase and quadrature bits
        det_inphase_bits=(np.real(sin_tap_equalisation)>0)
        det_quad_bits=(np.imag(sin_tap_equalisation)>0)

        # BER Simulation
        BER[j]=BER[j]+np.sum(det_inphase_bits!=inphase_bits)+np.sum(det_quad_bits!=quad_bits)

BER=BER/nSubcarriers/2/nblocks

# Theoretical BER calculation
# Effective SNR
effective_SNR = (L*SNR)/nSubcarriers



# BER Theoretical
BER_theory = 0.5*(1-np.sqrt(effective_SNR/(2+effective_SNR)))


# Superimposed plot BER theoretical and BER simulation
plt.figure(1)
plt.yscale('Log')
plt.plot(SNR_db,BER,'r')
plt.plot(SNR_db,BER_theory,'k:o')
plt.grid(1,which='both')
plt.suptitle('BER for OFDM wireless system')
plt.xlabel('SNR(in dbs)')
plt.ylabel('BER')
plt.legend(["BER Simulation","BER THEORY"])
plt.show()



