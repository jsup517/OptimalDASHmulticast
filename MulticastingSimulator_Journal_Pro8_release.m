%% Video Multicasting over LTE
clear all; close all; clc;
%% Videos - SVC or DASH
Rv = [1.071e6 1.662e6 2.617e6 3.305e6;...
      1.055e6 1.568e6 2.186e6 3.127e6;...
      1.011e6 2.828e6 3.676e6 4.376e6;...
      1.038e6 1.371e6 2.182e6 3.679e6];
PSNR = [35.1631 37.3271 38.2429 40.4971;...
      35.1631 37.3271 38.2429 40.4971;...
      35.1631 39.5843 41.7745 42.5728;...
      35.1631 37.3271 38.2429 41.7745];
[V L] = size(Rv); % V - num of video, L - num of layer
A = 1e5; alpha = 1./log(max(Rv')/A); beta = max(Rv'/A);
for vv=1:V, u(vv,:) = alpha(vv)*log(beta(vv).*Rv(vv,:)./max(Rv(vv,:))); end % utility
for ll=1:L, if ll==1, uvl(:,ll) = u(:,ll); else uvl(:,ll) = u(:,ll)-u(:,ll-1); end, end % utility difference

%% LTE parameters
nRB = 55%[15 25 35 45 55];%[6 15 25 50 100]; % Number of resource block in one OFDMA symbol (total : 1000)

%% Simulation parameters
SNRdBth = [1.4 2.24 3.54 5.6 8.2 11.01 13.81 17.92 22.4 25.2 30.8 36.4 42 47.6 52.08]; % Threshold SNR for 15-MCS
c = 12*7*[0.1523 0.2344 0.3770 0.6016 0.8770 1.1758 1.4766 1.9141 2.4063 2.7305 3.3223 3.9023 4.5234 5.1152 5.5547]/(0.5e-3); % efficiency
c_MSML = 12*7*[0.1523 0.2344 0.3770 0.6016 0.8770 1.1758 1.4766 1.9141 2.4063 2.7305 3.3223 3.9023 4.5234 5.1152 5.5547]/(0.5e-3); % efficiency
SNRth = 10.^(SNRdBth/10);
N = 100; % FEC block size
fm = 0.91; % FEC coding margin
iter = 10; % iteration

FECblock = 100;
tc = 100;  lambda = 1.5;
%% Users
Dist = 10e-1;                 % cell size
MinBsMs = 36e-3;              % minimum BS-MS distance 
nDataSubCPerSlot = 12;

xlsSheet = 'fxt=3x2';
xlsFileName_Ped1024 = 'ITU Ped B-fft1024-dB new.xls';
xlsFileName_Ped2048 = 'ITU Ped B-fft2048-dB new.xls';
xlsFileName_Veh1024 = 'ITU Veh A-fft1024-dB new.xls';
xlsFileName_Veh2048 = 'ITU Veh A-fft2048-dB new.xls';
FDdatabase = xlsread(xlsFileName_Ped2048, xlsSheet);
FD = mean(FDdatabase,2);
% Distance generation
simNn = 400:20:400; 
simcount = 0;
for nRBN = nRB
    Ne = nRBN;
for simN = simNn
    simcount = simcount+1
    %% Proposed initialization
    Q = ones(simN*V,length(SNRdBth)); Qa = ones(simN*V,length(SNRdBth)); I = ones(simN*V, length(SNRdBth));
    %% Start packet transmission
    for sim=1:iter
    %% SNR generation - PRB (Rayleigh distribution)
        if mod(sim,10)==1
            for ui = 1:simN*V % SNR generation for the 100-OFDMA frames (1-FEC block)
                di(ui) = sqrt(MinBsMs*MinBsMs+rand*(Dist*Dist-MinBsMs*MinBsMs));
                for ee=1:Ne*FECblock*(iter+1)
                    blockdB(ui,ee) = DistanceToSNR_LTE(di(ui), FD(ceil(1000*rand),:), nDataSubCPerSlot); 
                    block(ui,ee) = 10^(blockdB(ui,ee)/10);
                end
                for ui2=1:FECblock*(iter+1)
                    user(ui,ui2) = block(ui,Ne*(ui2-1)+1);%-lambda*log(sum(exp(-block(Ne*(ui2-1)+1:Ne*ui2)/lambda))/Ne);
                    userdB(ui,ui2) = blockdB(ui,Ne*(ui2-1)+1);%10*log10(user(ui,ui2));
                end
                pd1 = fitdist(userdB(ui,:)','normal');
            end
            for jj=1:length(SNRth), Q(:,jj) = sum(user(:,1:sim*FECblock)>SNRth(jj),2)/FECblock; end % Estimated receiving rates for Proposed
        end
        for jj=1:length(SNRdBth), Qa(:,jj) = sum(user(:,sim*FECblock+1:(sim+1)*FECblock)>SNRth(jj),2)/FECblock; end % Actual Receiving rate w.r.t MCS

        %% Optimal solution - Proposed
        mu = min(mean(userdB(:,1:sim*FECblock)'));
        pd = fitdist((mean(userdB(:,1:sim*FECblock)')-mu+0.01)','weibull');
        [RB MCS KN] = DASHMULNM6_DASH(pd.b, pd.a, mu, pd1.sigma, Rv, simN,nRBN, A) % simple search
        KN = KN+(RB<1);
        if sum(sum(KN>=1))==16, KN = 0.6*ones(4,4); end
        % Grouping
        for vv=1:V, for ll=1:L, group(sim,vv,find(Q(simN*(vv-1)+1:simN*vv,MCS(vv,ll))>KN(vv,ll))) = ll; end, end
        for vv=1:V
            for ll=L:-1:1
                groupindex = find(group(sim,vv,:)==ll);
                if isempty(groupindex), Nm(vv,ll,sim) = 0; else 
                    Nm(vv,ll,sim) = sum(Qa(simN*(vv-1)+groupindex,MCS(vv,ll))>=KN(vv,ll));     
                    if ll>1, group(sim,vv,groupindex(find(Qa(simN*(vv-1)+groupindex,MCS(vv,ll))<KN(vv,ll))))=ll-1; end
                end
            end
        end
        Utemp(sim) = sum(sum(Nm(:,:,sim).*u,2));
        PSNRtemp_pro(sim) = sum(sum(Nm(:,:,sim).*PSNR,2))/(4*400);
        %% Resource allocation for other schemes
        for b=1:FECblock
            % Proposed
            for jj=1:length(SNRdBth) 
                I(:,jj) = user(:,sim*FECblock+b)>SNRth(jj);
                Q(:,jj) = (1-1/tc)*Q(:,jj)+1/tc*I(:,jj); 
            end 
        end
    end

    %% Measure Packet Errors, Utility, Throughput
    U(simcount) = mean(Utemp(2:end));
    PSNR_pro(simcount) = mean(PSNRtemp_pro(2:end));

    finalNm(simcount,:) = sum(mean(Nm,3));
    finalNmL(:,:,simcount) = mean(Nm,3);
    finalNmp(simcount,:) = finalNm(simcount,:)/(simN*V);
    finalNm_min(simcount,:) = finalNmp(simcount,:)-min(reshape(sum(Nm),[L iter 1])')/(simN*V);
    finalNm_max(simcount,:) = finalNmp(simcount,:)-max(reshape(sum(Nm),[L iter 1])')/(simN*V);
    finalMCS(simcount,:) = mean(MCS);
    finalRB(simcount,:) = sum(mean(RB,3));
    finalRBL(:,:,simcount) = mean(RB,3);
    finalRB_min(simcount) = min(sum(sum(RB))); finalRB_max(simcount) = max(sum(sum(RB)));
    finalKN(simcount,:) = mean(KN);
    AvgThroughput(simcount) = sum(sum(Rv.*mean(Nm,3)))/(simN*V);
    Var(simcount,:) = mean(sum(((group(1:end-1,:,:)-group(2:end,:,:))~=0)),3)/iter;%mean(var(group),3);
    
end
end
%% Plots
figure(1);
subplot(1,2,1); grid on; hold; plot(nRB,U./(simNn*V),'-o','linewidth',2);
subplot(1,2,2); grid on; hold; plot(nRB,AvgThroughput'./sum(finalRB(:,:),2)/1e3,'-o','linewidth',2);
figure(2); hold; plot(nRB,PSNR_pro); 
figure(3); hold; plot(nRB,AvgThroughput'./sum(finalRB(:,:),2)/1e3); 
figure(4); hold; plot(nRB,U./simNn./sum(finalRB(:,:),2)'); 