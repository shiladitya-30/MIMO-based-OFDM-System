close all; clear all; rng('shuffle'); 
SNRdB = [1:5:56];
%subs
Nsub = 64;      %number of subcarriers
Ncp = round(Nsub/10);  %length of CP
numBlocks = 100;      
L = 3;       %number of channel taps between each transmit and receive antenna
r = 2; t = 2;   %number of receive and transmit antennas
BER = zeros(size(SNRdB));  
SNR = zeros(size(SNRdB));
for ix1 = 1:numBlocks
    ix1
    bits = randi([0,1],[t,Nsub]);  %message bits
    h = 1/sqrt(2)*(randn(r,t,L) + j*randn(r,t,L));   
    Hfreq = fft(h,Nsub,3);        
    ChNoise = randn(r,Nsub) + j*randn(r,Nsub);
    for ix2 = 1:length(SNRdB)
        SNR(ix2) = 10^(SNRdB(ix2)/10);
        RxbitsWithoutCP = zeros(r,Nsub);
        Loadedbits = sqrt(SNR(ix2)) * (2*bits - 1);   %BPSK symbols from message bits
        Txsamples = ifft(Loadedbits, Nsub, 2);        %transmitting symbols
        TxsamplesCP = [Txsamples(:, Nsub-Ncp+1:Nsub),Txsamples];   %transmitting symbols with CP
        for tx=1:t
            TxsamplesCp_txi = TxsamplesCP(tx,:);
            Rxbits = [];
            for rx = 1:r
                h_rxi_txi = squeeze(h(rx,tx,:));
                Rxbits = [Rxbits; conv(h_rxi_txi, TxsamplesCp_txi)];   %received bits with CP
            end
            RxbitsWithoutCP = RxbitsWithoutCP + Rxbits(:,Ncp+1:Ncp+Nsub);  %received bits without CP
        end
        RxbitsWithoutCP = RxbitsWithoutCP + ChNoise;
        RxbitsFFT = fft(RxbitsWithoutCP, Nsub, 2);
        for nx = 1:Nsub
            Hsub = Hfreq(:,:,nx);    
            ProcessedBits = pinv(Hsub)*RxbitsFFT(:,nx);
            DecodedBits = (real(ProcessedBits)>=0);     %decoding with threshold=0
            BER(ix2) = BER(ix2) + sum(DecodedBits ~= bits(:,nx));  %total bits in error
        end
    end
end
eSNR = L*SNR/Nsub;   %effective SNR
BER = BER/(numBlocks*Nsub*t);   
 
semilogy(SNRdB,BER,'b  s','linewidth',2.0);
hold on;
semilogy(SNRdB,0.5*(1-sqrt(eSNR./(2+eSNR))),'r -. ','linewidth',2.0);
axis tight;
grid on;
legend('MIMO ZF','Theory Single Ant')
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR(dB)');
