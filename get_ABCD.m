function [A,B,C,D] = get_ABCD(Es,Fs)
N=length(Es)-1; %order of the filter 
%ABCD matrix
for k=1:N+1
    if  mod(k,2)~=0
        A(N+2-k)=1i*imag(Es(N+2-k)+Fs(N+2-k));
        B(N+2-k)=real(Es(N+2-k)+Fs(N+2-k));
        C(N+2-k)=real(Es(N+2-k)-Fs(N+2-k));
        D(N+2-k)=1i*imag(Es(N+2-k)-Fs(N+2-k));
    end
    if  mod(k,2)==0   
        A(N+2-k)=real(Es(N+2-k)+Fs(N+2-k));
        B(N+2-k)=1i*imag(Es(N+2-k)+Fs(N+2-k));
        C(N+2-k)=1i*imag(Es(N+2-k)-Fs(N+2-k));
        D(N+2-k)=real(Es(N+2-k)-Fs(N+2-k));        
    end    
end
end