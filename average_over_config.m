clear
tic
variables;
w = w_00+w_D-4*w_0:0.0001:w_00+w_D+4*w_0;
Ncoh=zeros(ncon,1);
allshifts=zeros(ncon,N);
Stot = zeros(1,length(w));
Atot = zeros(1,length(w));
if ncon>1
 disp('disorder=yes!')
end 
results=zeros(ncon,N+2);
for index=1:ncon
 spectra_site;
 cal_Ncoh;
 if PBC==1
 CofS=coherencefT(coeff,energy,n1,v1,FC,kb,T,N,v_max,w_0,dim,PBC,threshold);
 if index==1
  C_CT = CofS;
 else
  C_CT = C_CT + CofS;
 end
 end
 Stot = Stot + S;
 Atot = Atot + A;
 allshifts(index,:)=siteshift';
 if PBC==1
 results(index,:)=[index,Ncoh(index),CofS];
 end
 index
end
disp('S.No.,Ncoh,C(s)(length=N)');
if PBC==1
results
end

Stot = Stot/ncon;
Atot = Atot/ncon;
if PBC==1
C_CT= C_CT/ncon; %coherence function averaged over T and disorder
end
if max(Stot)>0
ss=max(Stot);
Stot = Stot/ss;
end
if max(Atot)>0
aa=max(Atot);
Atot = Atot/aa;
%end
%save('siteshifts.mat','allshifts');
fileID=fopen('siteshifts.out','w');
formatspec='%15.8e\n';
fprintf(fileID,formatspec,allshifts(:,1));
fclose(fileID);
if ncon>1
 meanshift = mean(allshifts(:,1))
 stdshift = std(allshifts(:,1))
 skewness_out = skewness(allshifts(:,1))
 kurtosis_out = kurtosis(allshifts(:,1))-3
end

%% added by JN
formatspec='%15.8e  %15.8e\n';
fileID=fopen('sitebasis_emi.out','w');
fprintf(fileID,formatspec,[wn;Stot]);
fclose(fileID);
fileID=fopen('sitebasis_abs.out','w');
fprintf(fileID,formatspec,[wn;Atot]);
fclose(fileID);
%% %%%%%%%%%%%%%%%
end

disp('Time in minutes');
toc/60

%nodis=load('abs_thf_nodis.out');
%ls_nodis=nodis(:,2);
%expt=load('gfit_thf_thinfilm.out'); %this is abs
%wn_expt=expt(:,1);
%ls_wn_expt=expt(:,5);
%[maxm,maxm_idx]=max(ls_wn_expt);
%[maxm,maxm_idx_sim]=max(Atot);
%shiftt=wn_expt(maxm_idx)-wn(maxm_idx_sim);
%for i=1:size(wn,2)
% wn(i)=wn(i)+shiftt;
%end
%plot(wn,Atot,'-b','LineWidth',2);
%hold on
%plot(wn,ls_nodis,'-m','LineWidth',2);
%hold on
%plot(wn_expt,ls_wn_expt,'-k','LineWidth',2);
%legend('Sim','Sim-nodis','expt')
%hold off



