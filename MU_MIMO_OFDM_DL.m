%>>>>>>>>>>>>>>>>>>>>>MU-MIMO-OFDM-Downlink>>>>>>>>>>>>>>>>>>>>>>
clear all;
nSubC = 32;  % Total Subcarriers
nCP = round(nSubC/4); 
nUEs = 4; %Total Users
Ns=1;
Nt = 64; % Total BS Antennas
No=1;
var_db = [0,-9.7,-19.2,-22.8,-30]; %Ped-A Channel
%var_db = [1,1,1,1,1]; %Rayleigh Fading Channel
var = 10.^(var_db/10);
Ch_taps=length(var_db); % Total Channel Taps
SNR_db = -20:5:30;SNR = 10.^(SNR_db/10);
uu_itrns=5;lrg_itrns=5;sml_itrns=5;

magic=rng;
for precoder=1:3
rng(magic);
%User Locations>>>>>>>>>>>>
for uu=1:uu_itrns
min_dist=30;max_dist=1000;
USER_location=randperm(max_dist-min_dist+1,nUEs)+min_dist;

%Large Scale Fading>>>>>>>>>>>
for it=1:lrg_itrns
path_exp=3.8;mu_lognrm=3;sigma_lognrm=10^(1/10);
mu_nr=log10((mu_lognrm^2)/(sqrt((mu_lognrm^2)+(sigma_lognrm^2))));
sigma_nr_sq=log10(1+((sigma_lognrm^2)/(mu_lognrm^2)));
sigma_shadow=10;noise_dbm=-94;sigma_noise=(10^(-94./10))*1e-3;
median_ch_gain=10^(-3.453);
beta_dem = USER_location.^path_exp;
Beta_vec=(median_ch_gain.*(10.^((sqrt(sigma_shadow).*randn(1,nUEs))/10)))./beta_dem;
beta_values=Beta_vec./sigma_noise;
beta_save{it,:,precoder}=beta_values;
Beta=diag(beta_values);

for jj=1:sml_itrns
    hpre=(randn(nUEs,Nt,Ch_taps)+1j.*randn(nUEs,Nt,Ch_taps)).*(reshape(sqrt(var/2),1,1,[]));
    chnl{it,jj,precoder}=hpre;
    D=repmat(sqrt(Beta),[1,1,size(hpre,3)]);
    h=pagemtimes(D,hpre);

    %CFR>>>>>>>>>>>>
    Hm=fft(h,nSubC,3);
    
    for kk=1:length(SNR_db)
        [precoder,uu,it,jj,kk]
        p_users=SNR(kk)/nUEs;
        %Precoder>>>>>>>>>>>
        %MRT>>>>>>>>>>>>>>
        P1=zeros(Nt,nUEs,nSubC);
        if precoder==1 
            for m=1:nSubC
                P1=permute(conj(Hm),[2,1,3]);
            end
        %ZF>>>>>>>>>>>>>
        elseif precoder==2 
            for m=1:nSubC
                P1(:,:,m)=(Hm(:,:,m)')/(Hm(:,:,m)*(Hm(:,:,m)'));
            end
        %MMSE>>>>>>>>>>>>
        elseif precoder==3
            for m=1:nSubC
                P1(:,:,m)=(Hm(:,:,m)')/((Hm(:,:,m)*(Hm(:,:,m)')+(No/SNR(kk))*eye(nUEs)));
            end
        end

        %Normalization>>>>>>>>>
        normP1=sqrt(sum(abs(P1).^2, 1));
        P = bsxfun(@rdivide, P1, normP1);

        for m_dash=1:nSubC
            for k_dash=1:nUEs
                
                %>>>Signal Power>>>>>>>>>>>>>>>>>>>>>>>>>>>
                sgnpwr = abs(Hm(k_dash,:,m_dash)*P(:,k_dash,m_dash))^2;
                num_term = sgnpwr;

                %>>>Interference Power II>>>>>>>>>>>>>>>>>>>>>>>
                MUI=zeros(1,nUEs);
                for k=1:nUEs
                    if(k~=k_dash)
                        MUI(k) = Hm(k_dash,:,m_dash)*P(:,k,m_dash);
                    end
                end
                
                dem_term =abs(sum(MUI))^2;
                SINR = (p_users*num_term)/(dem_term*p_users+1);
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
SE(precoder,:)=sum(SE_uu,1)/uu_itrns;
end

plot(SNR_db,SE(1,:),'LineWidth', 1.5, 'Marker', 'square', 'MarkerSize', 6 , 'MarkerFaceColor','g');
hold on;
plot(SNR_db,SE(2,:),'LineWidth', 1.5, 'Marker', 'square', 'MarkerSize', 6 , 'MarkerFaceColor','c');
plot(SNR_db,SE(3,:),'LineWidth', 1.5, 'Marker', 'square', 'MarkerSize', 6 , 'MarkerFaceColor','r');
legend('MRT','ZF','MMSE','Location','northwest');
xlabel("Total Transmit Power");
ylabel("Downlink Sum Rate (bits/s/Hz)");
title("Ped-A Channel");
%title("Rayleigh Fading Channel");
grid on;



