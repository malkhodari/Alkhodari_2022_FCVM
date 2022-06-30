function [sF1,sM1,mphi,MtoF,fetal_beats,inst_fHR] = Phase_maternal_to_fetal(Maternal_Ratio,sM1,sF1,frr)

MtoF1=[]; 
index=1;
nmt=sM1;

mphi=[]; MtoF=[];
fetal_beats=0; fix=1;
inst_fHR=[];
index=1;
w1=1; w2=0;
ix=1;
%****This Loop to calculate maternal phase and relative Phase
for i=1:Maternal_Ratio:length(sM1)-Maternal_Ratio
    
    interval_m=nmt(i):1:nmt(i+Maternal_Ratio);
    tm=sF1( sF1< interval_m(end) & sF1>= interval_m(1)); 
    fetal_beats(fix)=length(tm);  fix=fix+1;
    g(ix)=length(interval_m);
    ix=ix+1;
    w2=w2+length(tm);
    
    cnstp=length(nmt(1):1:nmt(1+Maternal_Ratio));
%     
%     ist(1:length(tm))=mean(frr(w1:w2)); w1=w1+length(tm);
%     inst_fHR=[inst_fHR;ist'];
%     ist=[];
M_beats_interval=nmt(i+Maternal_Ratio)-nmt(i);
 
%******This to calculate Maternal Phase**********************
t=interval_m(1:end-1);
m_tmp=(2.*pi.*(t-nmt(i))./M_beats_interval)+(2*pi*0);
mphi=[mphi;m_tmp'];
%************************************************************

%******This to calculate Relative Phase**********************
 RelativeMF=(2.*pi.*(tm-nmt(i))./M_beats_interval)+(2*pi*0); %ORIGINAL PHASE
%  RelativeMF=exp(0.5.*pi.*(tm-nmt(i))./M_beats_interval)+(2*pi*0);
MtoF=[MtoF;RelativeMF'];
%************************************************************
MtoF_a{index}=RelativeMF./(2*pi);
% MtoF1(index)=median(RelativeMF./(2*pi));
% MtoF2(index)=max(RelativeMF./(2*pi));
index=index+1;
end
%*************************End of Loop*******************************

%    MtoF=MtoF./(2*pi);
end

