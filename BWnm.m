function NRF = BWnm(uy,K,alpha,n,beta,gama,deltaA,deltaV,deltaNu,pw,qw,zw,Fi,deltaFi,landa,dq)
% This function calculated the hysteresis part of Bouc-wen given
% the reletive displacement and velocity at the connection point
% zdot=qdot-beta*abs(qdot)*(abs(z)^(n-1))*z-gama*qdot*(abs(z)^n);

% Copyright 2022, All Rights Reserved
% Coded by Sina Safari
% University of Exeter
% ---------------------------------------------------------------
global HISTbw

% Hysteretic effects
A=1-deltaA*HISTbw.E;
V=1+deltaV*HISTbw.E;
Nu=1+deltaNu*HISTbw.E;
%Force terms
HISTbw.kbw=(A - V.*(beta.*sign(HISTbw.z).*sign(dq) + gama).* abs(HISTbw.z).^(n))./Nu;
% dz = HISTbw.kbw.*dq./uy;
hz=(1-(zw*(1-exp(-pw*HISTbw.E)))*exp(-(HISTbw.z*sign(dq)-qw*((1/((1+deltaV*HISTbw.E)*(beta+gama)))^(1/n)))^2/((Fi+deltaFi*HISTbw.E)*(landa+(zw*(1-exp(-pw*HISTbw.E)))))^2));
dz = hz*(HISTbw.kbw.*dq./uy);
LinearTerm = alpha.*K.*dq;
HystereticTerm = (1-alpha).*K.*uy.*dz;
HISTbw.K = alpha.*K + (1-alpha).*K.*HISTbw.kbw;
HISTbw.E=HISTbw.E+HISTbw.z*dq;

HISTbw.z=HISTbw.z+dz;
HISTbw.NRF=HISTbw.NRF+LinearTerm+HystereticTerm-K.*dq;
NRF=HISTbw.NRF;