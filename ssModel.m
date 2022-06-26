function dqdt = ssModel(t,q,para,P,ts)

M=para.M;
C=para.C;
K=para.K;
beta=para.beta;
gama=para.gama;
alpha=para.alpha;
n=para.n;
uy=para.uy;

deltaA=para.deltaA;
deltaV=para.deltaV;
deltaNu=para.deltaNu;

% Bouc-Wen pinching parameters
pw=para.pw;
qw=para.qw;
zw=para.zw;
Fi=para.Fi;
deltaFi=para.deltaFi;
landa=para.landa;

P=interp1(ts,P,t);
% q=q(1);qdot=q(2);z=q(3),e=q(4);
dqdt = [q(2);
        M\(P-C*q(2)-alpha*K*q(1)-(1-alpha)*uy*K*q(3));
        (1-(zw*(1-exp(-pw*q(4))))*exp(-(q(3)*sign(q(2))-qw*((1/((1+deltaV*q(4))*(beta+gama)))^(1/n)))^2/((Fi+deltaFi*q(4))*(landa+(zw*(1-exp(-pw*q(4))))))^2))*...
        (1/uy)*(1/(1+deltaNu*q(4)))*((1-deltaA*q(4))*q(2)-(1+deltaV*q(4))*(beta*abs(q(2))*(abs(q(3))^(n-1))*q(3)+gama*q(2)*(abs(q(3))^n)));
        q(2).*q(3)];