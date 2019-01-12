function [CQ]=get_Quartet(M,N)
a=M(3,8)*M(4,5)*M(5,6)-M(3,4)*M(6,7)*M(7,8)+M(3,8)*M(4,7)*M(6,7);
b=M(3,4)*M(4,7)*M(7,8)-M(3,8)*(M(4,5)^2-M(5,6)^2-M(6,7)^2+M(4,7)^2);
c=-M(3,8)*(M(4,7)*M(6,7)+M(4,5)*M(5,6));
t1=roots([a b c]);
Tk=[4 5 3];
Tl=[7 8 6];
Tm=[4 5 5];
Tn=[5 6 6];
Tc=[-1 -1 1];
i=[4 5 6 3];
j=[6 7 8 5];
theta_t=atan(t1(1));
CQ=M;
R=eye(N+2);
R(i(1),i(1))=cos(theta_t);
R(j(1),j(1))=cos(theta_t);
R(j(1),i(1))=sin(theta_t);
R(i(1),j(1))=-sin(theta_t);
CQ=R*CQ*R';
for in=2:length(j)
    theta_t=Tc(in-1)*atan(M(Tk(in-1),Tl(in-1))/M(Tm(in-1),Tn(in-1)));
    R=eye(N+2);
    R(i(1),i(1))=cos(theta_t);
    R(j(1),j(1))=cos(theta_t);
    R(j(1),i(1))=sin(theta_t);
    R(i(1),j(1))=-sin(theta_t);
    CQ=R*CQ*R';
end
display('Cascaded Quartet');
printMatrix(CQ);
end