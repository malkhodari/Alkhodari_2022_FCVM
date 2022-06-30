function [lam]=lambda_A(psi,N)
    %24/10
    %Lambda is used to quantize the phase shift pattern of two ecg signals
    %use phase shift function to get psi
    %feed psi to lambda to get the lambda pattern
    %N is the size of the computation window, the higher it is the smoother the
    %output.
    p=1;lam=[];
    %psi=mod(psi,2*pi);
    while p<length(psi)-N
        %sum_= 0;
        q=[];sumc=0; sums=0;
        for q=p:p+N-1
            if psi(q)==0
                a=0;
            else
            a=exp(1i*psi(q));
            end
            sumc=sumc+real(a);
            sums=sums+imag(a);
            %sum_=sum_+exp(1i*psi(q));
        end
        lam(p)=(sumc/N).^2 + (sums/N).^2;
        %lam(p)=(real(sum_)/N)^2 ;      
        p=p+1;
    end
    lam=medfilt1((lam),5);
end