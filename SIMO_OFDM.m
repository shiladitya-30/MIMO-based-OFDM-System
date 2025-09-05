clear all; close all; clc;
SNR_dB = [1:5:45];
SNR = 10.^(SNR_dB/10);
BER = zeros(length(SNR),1);
 
Nsub = 64;   %Number of subcarriers
Ncp = round(Nsub/10);  %Number of samples added as CP
 
numTaps = 3;   %Number of channel taps
numBlocks = 10000;
numAnt = 2;    %Number of receive antennas
 
for L = 1:numBlocks
    Bits = randi([0,1],1,Nsub);    %Generating message bits
    h = 1/sqrt(2) * (randn(numAnt,numTaps) + 1j*randn(numAnt,numTaps));  
    Hfreq = fft(h, Nsub, 2);  %columnwise N point FFT of channel taps
    ChNoise = randn(numAnt, Nsub+Ncp+numTaps-1) + 1j*randn(numAnt, Nsub+Ncp+numTaps-1);
    for K = 1:length(SNR)
        LoadedBits = sqrt(SNR(K))*(2*Bits - 1);  %BPSK symbols
        TxSamples = ifft(LoadedBits);    %IFFT of symbols at transmitter
        TxSamples_CP = [TxSamples(Nsub-Ncp+1:Nsub) TxSamples];   %Adding CP 
        
        Rxbits = [];
        for rxi = 1:numAnt
            Rxbits = [Rxbits; conv(TxSamples_CP, h(rxi,:)) + ChNoise(rxi,:)];
        end
        
        RxbitsWithoutCP = Rxbits(:, Ncp+1:Ncp+Nsub);    %Removal of CP
        RxbitsFFT = fft(RxbitsWithoutCP, Nsub, 2);  %FFT of received symbols at receiver
        ProcessedBits = sum(RxbitsFFT.*conj(Hfreq));  %Matched filtering at receiver
        DecodedBits = (real(ProcessedBits)>=0);    %Decoding received bits
        BER(K) = BER(K) + sum(DecodedBits ~= Bits);  %BER for Kth SNR value
        
    end
end
 
BER = BER/(numBlocks*Nsub);
semilogy(SNR_dB, BER, 'b - s', 'linewidth', 2.0);
hold on; axis tight; grid on;
 
SNR_eff = numTaps/Nsub * SNR;
semilogy(SNR_dB,nchoosek(2*numAnt-1,numAnt-1)*1/2^numAnt./SNR_eff.^numAnt,'g -. ','linewidth',2.0);
 
xlabel('SNRdB')
ylabel('BER')
title('BER vs SNR for SIMO OFDM System');
legend('Simulation', 'Theory')
