clear all; 
close all; 
clc;
SNR_dB = [1:5:60];  %SNR range in dB
SNR = 10.^(SNR_dB/10);  %SNR in linear scale
BER = zeros(length(SNR),1);
 
Nsub = 64;  %Number of subcarriers
Ncp = round(Nsub/10);   %Number of samples added as CP
 
numTaps = 3;    %Number of channel taps
numBlocks = 10000;
 
for L = 1:numBlocks 
    Bits = randi([0,1],1,Nsub);    %Generating message bits
    h = 1/sqrt(2) * (randn(1,numTaps) + 1j*randn(1,numTaps));   
    Hfreq = fft(h, Nsub);    %N point FFT of Channel taps
    ChNoise=randn(1, Nsub+Ncp+numTaps-1) + 1j*randn(1, Nsub+Ncp+numTaps-1);

    for K = 1:length(SNR)
        LoadedBits = sqrt(SNR(K))*(2*Bits - 1);    %BPSK symbols
        TxSamples = ifft(LoadedBits);       %IFFT of symbols at transmitter
        TxSamples_CP = [TxSamples(Nsub-Ncp+1:Nsub) TxSamples];   %Adding CP 
        Rxbits = conv(TxSamples_CP, h) + ChNoise;    
        RxbitsWithoutCP = Rxbits(Ncp+1:Ncp+Nsub);   %Removal of CP
        RxbitsFFT = fft(RxbitsWithoutCP);     %FFT of received symbols at receiver
        ProcessedBits = RxbitsFFT./Hfreq;     %Equalization at receiver
        DecodedBits = (real(ProcessedBits)>=0);
        BER(K) = BER(K) + sum(DecodedBits ~= Bits);   %BER for Kth SNR value
        
    end
end
 
BER = BER/(numBlocks*Nsub);
semilogy(SNR_dB, BER, 'b - s', 'linewidth', 2.0);
hold on; 
axis tight; 
grid on;
 
SNR_eff = numTaps/Nsub * SNR;
semilogy(SNR_dB, (1/2*(1-sqrt(SNR_eff./(SNR_eff+2)))), 'r  s', 'linewidth', 2.0);
xlabel('SNRdB')
ylabel('BER')
title('BER vs SNR for SISO OFDM System');
legend('Simulation', 'Theory')
