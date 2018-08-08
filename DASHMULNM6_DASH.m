function [RBout MCSout KNout] = DASHMULNM6_DASH(alpha, beta, mu, sigma, Rv, Nt, N_RB, A),
%% DASH Multicast 
N = 100; % Num of Frame : 1 Frame = 10msec (total 1sec)
fm = 0.91; % Decoding margin
K = 1:N*fm; % Num of Data
SNRdBth = [1.4 2.24 3.54 5.6 8.2 11.01 13.81 17.92 22.4 25.2 30.8 36.4 42 47.6 52.08];
SNRth = 10.^(SNRdBth/10);
c = 12*7*[0.1523 0.2344 0.3770 0.6016 0.8770 1.1758 1.4766 1.9141 2.4063 2.7305 3.3223 3.9023 4.5234 5.1152 5.5547]/(0.5e-3);
nRB = 0:N_RB; % Num of Resource Blocks
nv = length(Rv); [Vn Nv] = size(Rv);
tempA = max(Rv')'; %A = 1e4; 
alpha1 = 1./log(tempA/A); beta1 = tempA/A;
for uu=1:length(alpha1), oriu(uu,:) = alpha1(uu).*log(beta1(uu).*Rv(uu,:)./tempA(uu)); end
for ll=1:Nv, 
    if ll==1, uvl(:,ll) = oriu(:,ll); 
    else uvl(:,ll) = oriu(:,ll)-oriu(:,ll-1); end 
end % utility difference

testRB = round(N_RB/Vn);
for pp = 1:Nv,
    for vv = 1:Vn,
        for ii = 1:length(c),
            KN = Rv(vv,pp)./(nRB*c(ii)); KN=KN.*(KN<1)+(KN>=1);% Minimum coding rate
            Am = 2*Rv(vv,pp)./(fm*c(ii));
            temp1 = 1-Am./nRB; temp1 = temp1.*(temp1>-1)-(temp1<=-1);
            temp2 = ((SNRdBth(ii)-mu-sigma*sqrt(2)*erfinv(temp1))./beta).^alpha;
            Nm{vv,pp}(ii,:) = Nt*exp(-temp2.*(temp2>0));
            dNm{vv,pp}(ii,1) = 0; dNm{vv,pp}(ii,2:length(Nm{vv,pp}(ii,:))) = Nm{vv,pp}(ii,2:end)-Nm{vv,pp}(ii,1:end-1);
        end
        DNm{vv,pp} = Nm{vv,pp}.*dNm{vv,pp};
        [temval temindex]=max(max(Nm{vv,pp}));
        [mval bMCS(vv,pp)]=max(Nm{vv,pp}(:,temindex));
        testNm(vv,pp) = (sum(Nm{vv,pp}(:,testRB)>(Nt-1)))>0;
    end
    if sum(testNm(:,pp)) == Vn, startNv(pp) = 1; end
end
if false
    comb(1,1:Vn) = Nv*ones(1,Vn);
    optRB = zeros(1,Vn);
    optKN = zeros(1,Vn);
    optMCS = zeros(1,Vn);
    optRv = zeros(1,Vn);
    optu = zeros(1,Vn);
else
    for uu=1:(Nv-1+1)^Vn
        for vv=1:Vn
            comb(uu,vv) = mod(round(uu/(Nv-1+1)^(vv-1)),(Nv-1+1))+1;
        end
    end
    optRB = zeros((Nv-1)^Vn,(Nv-1)*Vn);
    optKN = zeros((Nv-1)^Vn,(Nv-1)*Vn);
    optMCS = zeros((Nv-1)^Vn,(Nv-1)*Vn);
    optRv = zeros((Nv-1)^Vn,(Nv-1)*Vn);
    optu = zeros((Nv-1)^Vn,(Nv-1)*Vn);
end
[acomb bcomb] = sort(sum(comb,2));
comb = comb(bcomb,:);
t_pro = cputime;
for ll=1:1
    clear newRv; clear MCS; clear comv; clear u; clear U; clear optNm; clear combi; clear newIndex; clear MCS; clear tempMCS; combsize = 1;
    for vv=1:Vn
        combRv{vv} = nchoosek(Rv(vv,:),comb(ll,vv));
        combIndex{vv} = nchoosek(1:Nv,comb(ll,vv));
        tempu1 = nchoosek(oriu(vv,:),comb(ll,vv));
        tempu2 = zeros(size(tempu1,1),size(tempu1,2));
        tempu2(:,2:end) = tempu1(:,1:end-1);
        combu{vv} = tempu1-tempu2;%uvl(vv,startIndex:comb(ll,vv));
        combcount(vv) = size(combRv{vv},1);
        combsize = combsize * combcount(vv);
    end
    combi = makeComb(combcount);
    for kk=1:size(combi,1)
        vcount = 0; vstart = 1;
        for vv=1:Vn
            vcount = vcount + comb(ll,vv);
            newRv(kk,vstart:vcount) = combRv{vv}(combi(kk,vv),:);%nchoosek(Rv(vv,:),comb(ll,vv));
            newIndex(kk,vstart:vcount) = combIndex{vv}(combi(kk,vv),:);
            u(kk,vstart:vcount) = combu{vv}(combi(kk,vv),:);%tempu1-tempu2;%uvl(vv,startIndex:comb(ll,vv));
            comv(kk,vstart:vcount) = vv;
            vstart = vcount+1;
        end
        nRBm = round(N_RB/length(newRv(kk,:)))*ones(1,length(newRv(kk,:)));
        nRBm = floor(N_RB*[1/2 1/4 1/4 1/4]);
        nRBm = nRBm + (nRBm==0);
        mm=1; stopMCS = 1; MCS(kk,:) = ones(1,length(newRv(kk,:)));
        while stopMCS
            for rr=1:length(newRv(kk,:))
                [mval tempMCS(rr)]=max(Nm{comv(kk,rr),newIndex(kk,rr)}(:,nRBm(rr)));
                if tempMCS(rr)<bMCS(comv(kk,rr),newIndex(kk,rr)), tempMCS(rr)=bMCS(comv(kk,rr),newIndex(kk,rr)); end
            end
            epsilon = 1;
            error = 100; %nRBm = zeros(1,length(newRv(kk,:)));
            temp_min = -5; temp_max = 5;
            early_exit = 0; merror = 10100;
            while epsilon<error,
                early_exit = early_exit + 1;
                temp_mean = (temp_max+temp_min)/2;
                lcount = 1-1;
                for rr=1:length(newRv(kk,:)),
                    lcount = lcount + 1;
                    [min_val(rr) nRBm(rr)]=min(abs(temp_mean-log(u(kk,rr)*dNm{comv(kk,rr),newIndex(kk,rr)}(tempMCS(rr),:))));
                    if lcount == comb(ll,comv(kk,rr)), lcount = 1-1; end
                end
                if N_RB-sum(nRBm)>0, temp_max = temp_mean; else temp_min = temp_mean; end
                error = abs(N_RB-sum(nRBm));
                if error<merror, bRBm = nRBm; merror = error; end
                if early_exit>10, break; end
            end    
            if sum(bRBm)>N_RB, bRBm = ceil(N_RB*bRBm/sum(bRBm)); end
            lcount = 1-1;
            for rr=1:length(newRv(kk,:))
                lcount = lcount + 1;
                optNm(kk,rr) = Nm{comv(kk,rr),newIndex(kk,rr)}(tempMCS(rr),bRBm(rr));
                if lcount == comb(ll,comv(rr)), lcount = 1-1; end
            end
            if mm>2 
                if MCS(kk,:)==tempMCS, stopMCS = 0; 
                else MCS(kk,:) = tempMCS; end
            end
            mm=mm+1;
            if mm>10, break; end
        end
        opt_nRB(kk,1:length(bRBm)) = bRBm;
        U(kk) = 0;
        for pp=1:length(newRv(kk,:))
            U(kk) = U(kk) + u(kk,pp)*(optNm(kk,pp));   
        end
        U(kk) = U(kk).*(U(kk)<10000);
    end
    [optU(ll) opti(ll)] = max(U);
    optRB(ll,1:length(newRv(kk,:))) = opt_nRB(opti(ll),1:length(newRv(kk,:)));
    optMCS(ll,1:length(newRv(kk,:))) = MCS(opti(ll),1:length(newRv(kk,:)));
    optKN(ll,1:length(newRv(kk,:))) = newRv(opti(ll),:)./(opt_nRB(opti(ll),1:length(newRv(kk,:))).*c(MCS(opti(ll),1:length(newRv(kk,:)))));
    optRv(ll,1:length(newRv(kk,:))) = newRv(opti(ll),:);
    optu(ll,1:length(newRv(kk,:))) = u(opti(ll),:);
    optv(ll,1:length(newRv(kk,:))) = comv(opti(ll),:);
    optIndex(ll,1:length(newRv(kk,:))) = newIndex(opti(ll),:);
    optN(ll,1:length(newRv(kk,:))) = optNm(opti(ll),1:length(newRv(kk,:))); 
    optRB(ll,:);
end
e_pro = cputime-t_pro;

perfectU = find(optU>=Vn*Nt);
if isempty(perfectU)
    [maxU maxi] = max(optU);
else
    [minRB minRBi] = min(sum(optRB(perfectU,:),2));
    maxi = perfectU(minRBi);
end
RBout = zeros(Vn, Nv); MCSout = length(c)*ones(Vn, Nv); KNout = zeros(Vn, Nv);
lcount = 1-1;
for ii=1:sum(optRB(maxi,:)>0) 
    lcount = lcount + 1;
    RBout(optv(maxi,ii),optIndex(maxi,ii)) = optRB(maxi,ii);
    MCSout(optv(maxi,ii),optIndex(maxi,ii)) = optMCS(maxi,ii);
    if optKN(maxi,ii)==0, KNout(optiv(maxi,ii),optIndex(maxi,ii)) = 1;
    else KNout(optv(maxi,ii),optIndex(maxi,ii)) = optKN(maxi,ii); end
    if lcount == comb(maxi,optv(maxi,ii)), lcount = 1-1; end
end
figure(1); hold on;
lcount = 1-1;
for ii=1:sum(optRB(maxi,:)>0)
    lcount = lcount + 1;
    plot(optu(maxi,ii)*dNm{optv(maxi,ii),optIndex(maxi,ii)}(optMCS(maxi,ii),:));
    plot(optRB(maxi,ii),optu(maxi,ii)*dNm{optv(maxi,ii),optIndex(maxi,ii)}(optMCS(maxi,ii),optRB(maxi,ii)),'o'); 
    if lcount == comb(maxi,optv(maxi,ii)), lcount = 1-1; end
end
figure(2); hold on;
lcount = 1-1;
for ii=1:sum(optRB(maxi,:)>0)
    lcount = lcount + 1;
    plot(Nm{optv(maxi,ii),optIndex(maxi,ii)}(optMCS(maxi,ii),:));
    plot(optRB(maxi,ii),Nm{optv(maxi,ii),optIndex(maxi,ii)}(optMCS(maxi,ii),optRB(maxi,ii)),'o'); 
    if lcount == comb(maxi,optv(maxi,ii)), lcount = 1-1; end
end
[a b] = sort(sum(comb,2));
%figure(3); plot(a,optU(b),'+');