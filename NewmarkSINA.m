function [qdd,qd,q,z,E]=NewmarkSINA(Time,M,C,K,P,Idisp,Ivel,gama,beta,tol,maxit,p,NLF,rho,LL,T,FI)
% Copyright 2022, All Rights Reserved
% Coded by Sina Safari
% University of Exeter
% ---------------------------------------------------------------
clear global HISTbw
global HISTbw

fnorm=norm(P);
dt=Time(2)-Time(1);
i=0;
B=FI'*LL'*T'*rho';
q(:,i+1)=Idisp;
qd(:,i+1)=Ivel;
HISTbw.z=q(:,i+1);
HISTbw.E=HISTbw.z*qd(:,i+1).*dt;
HISTbw.NRF=zeros(size(q(:,i+1)));
NRF(:,1)=NLF(p,B,1,q(:,i+1),qd(:,i+1),qd(:,i+1),Idisp);
z(:,i+1)=HISTbw.z;
E(:,i+1)=HISTbw.E;
RF(:,1)=C*qd(:,1)+K*q(:,1)+NRF(:,1); 
qdd(:,1)=M\(P(:,1)-RF(:,1));
r(:,1)=P(:,1)-M*qdd(:,1)-RF(:,1);
k_eff=K+(gama/(beta*dt))*C+(1/(beta*(dt^2)))*M;

% Time incrementation  
for i=1:length(Time)-1
f_eff(:,i)=(P(:,i+1)-P(:,i))+((1/(2*beta))*M+((gama*dt)/(2*beta))*C-dt*C)*qdd(:,i)+((1/(beta*dt))*M+(gama/beta)*C)*qd(:,i);
Dq(:,i)=k_eff\f_eff(:,i);
    q(:,i+1)=q(:,i)+Dq(:,i);
    qd(:,i+1)=(gama/(beta*dt))*Dq(:,i)+(1-(gama/beta))*qd(:,i)+dt*(1-(gama/(2*beta)))*qdd(:,i);
        HISTbw.z=z(:,i);
        HISTbw.E=E(:,i);
    NRF(:,i+1)=NLF(p,B,1,q(:,i+1),qd(:,i+1),qd(:,i+1),qd(:,i+1)*dt);
    RF(:,i+1)=C*qd(:,i+1)+K*q(:,i+1)+NRF(:,i+1);
    qdd(:,i+1)=M\(P(:,i+1)-RF(:,i+1));
    f_eff(:,i)=f_eff(:,i)-NRF(:,i+1)+NRF(:,i);
    Dq(:,i)=k_eff\f_eff(:,i);
    r(:,i+1)=P(:,i+1)-M*qdd(:,i+1)-RF(:,i+1);

j=1;
while    norm(r(:,i+1))>tol&&(j<=maxit)
    q(:,i+1)=q(:,i)+Dq(:,i);
    qd(:,i+1)=(gama/(beta*dt))*Dq(:,i)+(1-(gama/beta))*qd(:,i)+dt*(1-(gama/(2*beta)))*qdd(:,i);
        HISTbw.z=z(:,i);
        HISTbw.E=E(:,i);
    NRF(:,i+1)=NLF(p,B,1,q(:,i+1),qd(:,i+1),qd(:,i+1),qd(:,i+1)*dt);
    RF(:,i+1)=C*qd(:,i+1)+K*q(:,i+1)+NRF(:,i+1);
    qdd(:,i+1)=M\(P(:,i+1)-RF(:,i+1));
    r(:,i+1)=P(:,i+1)-M*qdd(:,i+1)-RF(:,i+1);
    f_eff(:,i)=f_eff(:,i)-NRF(:,i+1)+NRF(:,i);
    Dq(:,i)=k_eff\f_eff(:,i);
    
    j=j+1;
   
end
z(:,i+1)=HISTbw.z;
E(:,i+1)=HISTbw.E;
end
assignin('base','HISTbw',HISTbw)