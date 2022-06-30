function [phiall, phi, tim]=phase_shift_A(F_ECG,M_ECG,fb,mb,t)
    %mb=number of maternal beats
    %fb=number of fetal beats
    %F_ECG= R peak locations of FECG in sec
    %M_ECG= R peak locations of MECG in sec
    % t is the time (sec) of MECG signals

    M_ECG=single(M_ECG*1000);
    F_ECG=single(F_ECG*1000);

    phiall=[];phi=[];tim=[]; phitemp=[];phitempall=[];ti=[];

    for k=1:fb:length(F_ECG)-fb
        temp=M_ECG(M_ECG>=F_ECG(k) & M_ECG <F_ECG(k+fb));
        ti=t(temp);

        if length(temp)==mb
            phitemp=2*pi*(F_ECG(k)-temp)/(F_ECG(k)-F_ECG(k+fb))+2*pi*k;
        else
           phitemp=zeros(1,length(temp));
        end

        phitempall=2*pi*(F_ECG(k)-temp)/(F_ECG(k)-F_ECG(k+fb))+2*pi*k;
        phiall=[phiall phitempall];
        phi=[phi phitemp];

        tim = [tim ti];
    end
    phi=mod(phi,2*pi)/(2*pi);
    phiall=mod(phiall,2*pi)/(2*pi);
end