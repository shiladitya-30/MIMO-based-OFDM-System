%>>>>>>>>>>>>>>>>>>>>>Massive-MIMO-OFDM-Uplink>>>>>>>>>>>>>>>>>>>>>>
clear all;
load loc_save.mat
nSubC = 32;  %Total Subcarriers
nCP = round(nSubC/4);Ns=1;
nUEs = 4; %Total Users
nNrs = [64,256]; %Total BS Antennas
Nr=nNrs(2);
No=1;
var_db = [0,-9.7,-19.2,-22.8,-30]; %Ped-A Channel
%var_db = [1,1,1,1,1]; %Rayleigh Fading Channel
var = 10.^(var_db/10);
Ch_taps=length(var_db); %Total Channel Taps
SNR_db = -20:5:30;SNR = 10.^(SNR_db/10);
lrg_itrns=5;sml_itrns=5;uu_itrns=5;

for it=1:lrg_itrns
    for jj=1:sml_itrns
        htotl{it,jj}=(randn(nNrs(2),nUEs,Ch_taps)+1j.*randn(nNrs(2),nUEs,Ch_taps)).*(reshape(sqrt(var/2),1,1,[]));
    end
end

magic=rng;
for combiner=1:3
rng(magic);
%User Locations>>>>>>>>>>>>
for uu=1:uu_itrns
min_dist=30;max_dist=1000;
USER_location=loc_save{uu,:};
%Large Scale Fading>>>>>>>>>>>>>>>
for it=1:lrg_itrns
    path_exp=3.8;mu_lognrm=3;sigma_lognrm=10^(1/10);
    mu_nr=log10((mu_lognrm^2)/(sqrt((mu_lognrm^2)+(sigma_lognrm^2))));
    sigma_nr_sq=log10(1+((sigma_lognrm^2)/(mu_lognrm^2)));
    sigma_shadow=10;noise_dbm=-94;sigma_noise=(10^(-94./10))*1e-3;
    median_ch_gain=10^(-3.453);
    beta_dem = USER_location.^path_exp;
    Beta_vec=(median_ch_gain.*(10.^((sqrt(sigma_shadow).*randn(1,nUEs))/10)))./beta_dem;
    beta_values=Beta_vec./sigma_noise;
    beta_save{it,:}=beta_values;
    Beta=diag(beta_values);

    %Small Scale Fading>>>>>>>>>>>>>>
    for jj=1:sml_itrns
        htemp=htotl{it,jj};
        hpre=htemp(1:Nr,:,:);
        D=repmat(sqrt(Beta),[1,1,size(hpre,3)]);
        h=pagemtimes(hpre,D);
        %CFR>>>>>>>>>>>>
        Hm=fft(h,nSubC,3);

        for kk=1:length(SNR_db)
            [combiner,uu,it,jj,kk]
            %User Power Allocation>>>>>>>>>
            p_user=SNR(kk)/nUEs;%Equal Power
            % Combiner>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            W1=zeros(nUEs,Nr,nSubC);
            %MRC>>>>>>>
            if combiner==1  
                for m=1:nSubC
                    W1(:,:,m)=Hm(:,:,m)';   
                end
            %ZF>>>>>>>
            elseif combiner==2  
                for m=1:nSubC
                    W1(:,:,m)=inv((Hm(:,:,m)')*Hm(:,:,m))*(Hm(:,:,m)');
                end
            %MSSE>>>>>>>>>
            elseif combiner==3 
               for m=1:nSubC
                    W1(:,:,m)=inv(((Hm(:,:,m)')*Hm(:,:,m))+((No/SNR(kk))*eye(nUEs)))*(Hm(:,:,m)');
               end
            end

            %Normalization>>>>>>>>>
            normW1=sqrt(sum(abs(W1).^2, 2));
            W = bsxfun(@rdivide, W1, normW1);
            
            for m_dash=1:nSubC
                for k_dash=1:nUEs
                    %Signal Power>>>>>>>>>>>>>>>>
                    sgnpwr=abs(W(k_dash,:,m_dash)*Hm(:,k_dash,m_dash))^2;
                    num_term=sgnpwr;
    
                    %Interference Power>>>>>>>>>>>>>>
                    MUI=zeros(1,nUEs);
                    for k=1:nUEs
                        if(k~=k_dash)
                            MUI(k)=W(k_dash,:,m_dash)*Hm(:,k,m_dash);
                        end
                    end
                    dem_term=abs(sum(MUI))^2;
    
                    SINR = (p_user*num_term)/(dem_term*p_user+1);
                    SUM_rate(k_dash) = log2(1+SINR);
                end
                SUM_rate_avr(m_dash)=sum(SUM_rate);
            end
        Srate(jj,kk) = sum(SUM_rate_avr)/nSubC;
        end
    end
SE_avr(it,:) = sum(Srate,1)/sml_itrns;
end
SE_uu(uu,:)=sum(SE_avr,1)/lrg_itrns;
end
SE(combiner,:)=sum(SE_uu,1)/uu_itrns;
end

plot(SNR_db,SE(1,:),'LineWidth', 1.5, 'Marker', 'square', 'MarkerSize', 6 , 'MarkerFaceColor','g');
hold on;
plot(SNR_db,SE(2,:),'LineWidth', 1.5, 'Marker', 'square', 'MarkerSize', 6 , 'MarkerFaceColor','b');
plot(SNR_db,SE(3,:),'LineWidth', 1.5, 'Marker', 'square', 'MarkerSize', 6 , 'MarkerFaceColor','c');
%legend('MRC','ZF','MMSE','Location','northwest');
xlabel("Total Transmit Power (dB)");
ylabel("Uplink Sum Rate (bits/s/Hz)");
title("Ped-A Channel");
%title("Rayleigh Fading Channel");
grid on;