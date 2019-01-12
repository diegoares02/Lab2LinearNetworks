clear;clc;
epsilon=43.4176;
RL=20;
Ps=[1i 0 4.85*1i 0 5.6644*1i];
Fs=[1 0 2.1412 0 1.4691 0 0.3369 0 0.0132];
Es=[1 1.9462 4.035 4.718 4.8027 3.4133 1.8758 0.6764 0.1322];
N=length(Es)-1;
[A,B,C,D]=get_ABCD(Es,Fs);
[M]=get_Matrix(A,B,C,D,Ps,epsilon,N);
[F]=get_Folded(M,N);
[CQ]=get_Quartet(F,N);
W=eye(N+2);
W(1,1)=0;
W(N+2,N+2)=0;
R2=zeros(N+2);
R2(1,1)=1;
R2(N+2,N+2)=1;
indice=1;
for w=-5:0.01:5
    Am1=-1i*R2+w*W+M;
    Am2=-1i*R2+w*W+F;
    Am3=-1i*R2+w*W+CQ;
    Am1=inv(Am1);
    Am2=inv(Am2);
    Am3=inv(Am3);
    S21t(indice)=2*1i*(Am1(N+2,1));
    S11t(indice)=-(1+2*1i*(Am1(1,1)));
    S21f(indice)=2*1i*(Am2(N+2,1));
    S11f(indice)=-(1+2*1i*(Am2(1,1)));
    S21q(indice)=2*1i*(Am3(N+2,1));
    S11q(indice)=-(1+2*1i*Am3(1,1));
    indice=indice+1;
end
d=-5:0.01:5;
S21tran=20*log10(abs(S21t));
S11tran=20*log10(abs(S11t));
S21fol=20*log10(abs(S21f));
S11fol=20*log10(abs(S11f));
S21qua=20*log10(abs(S21q));
S11qua=20*log10(abs(S11q));

figure(1)
plot(d,S21tran);
hold on
plot(d,S11tran,'r');
axis([min(d) max(d) -350 10])
title('Transversal Configuration')
xlabel('Frequency');
ylabel('Level(dB)');
grid on

figure(2)
plot(d,S21fol);
hold on
plot(d,S11fol,'r');
axis([min(d) max(d) -350 10])
title('Folded Configuration')
xlabel('Frequency');
ylabel('Level(dB)');
grid on

figure(3)
plot(d,S21qua);
hold on
plot(d,S11qua,'r');
axis([min(d) max(d) -350 10])
title('Quarter Configuration')
xlabel('Frequency');
ylabel('Level(dB)');
grid on