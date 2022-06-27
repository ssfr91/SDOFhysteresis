%% This Code is a Demo for validation of newmark solver dealing with
% hysteretic nonlinear models (Bouc-Wen-Barber-Noori)

% Copyright 2022, All Rights Reserved
% Coded by Sina Safari
% University of Exeter
% ---------------------------------------------------------------

clc
clear
close all
% Linear Parameters
BWpara.M=1;
BWpara.C=1;
BWpara.K=1000;
% Bouc-Wen main Parameters
BWpara.beta=2.5;
BWpara.gama=1;
BWpara.alpha=0.3;
BWpara.n=3;
BWpara.uy=1.5e-5;      % yielding displacement
% Bouc-Wen degradation parameters
BWpara.deltaA=0;     %
BWpara.deltaV=1e4;   %strength
BWpara.deltaNu=1e4;  %stiffness
% Bouc-Wen pinching parameters
BWpara.pw=1e5;
BWpara.qw=0.2;
BWpara.zw=0.9;
BWpara.Fi=0.2;
BWpara.deltaFi=1e4;
BWpara.landa=0.1;

% dynamic characteristics
wn=sqrt((BWpara.K)/BWpara.M);
fn=wn/2/pi;
dr=BWpara.C/(2*BWpara.M*wn);
%% Excitation
Fs=2^10;
dt=1/Fs;
ts=0:dt:10;

% % Symmetric Pulse
% t0=2;
% Tp=0.5;
% g=9.806; %m/s^2
% ap=0.8*g;
% u=para.M*ap*(1-(((2*pi^2)*(ts-t0).^2)/Tp^2)).*exp(-(1)*(pi^2)*((ts-t0).^2)/(Tp^2));

% % Sine Excitation
A = 1e-2;                 % excitation amplitude.
w=(2*pi)*1;
u=A*sin(w*ts);

%% Simulation (ode)
ndof=1;
options = odeset('AbsTol',[1e-10*ones(1,4*ndof)]);
IC=zeros(4*ndof,1); % Initial conditions
IC(1)=2e-5;
kbw=(1 - (BWpara.beta.*sign(0).*sign(IC(1)) + BWpara.gama).* abs(0).^(BWpara.n));
IC(3) = kbw.*IC(1)./BWpara.uy;
disp('ODE45')
tic
[Time,x] = ode45(@(t,q) ssModel(t,q,BWpara,u,ts),ts,IC,options);
toc
Time=ts;
q=x(:,1)';qdot=x(:,2)';z=x(:,3)';e=x(:,4)';
RF=BWpara.alpha*BWpara.K*q+(1-BWpara.alpha)*BWpara.uy*BWpara.K*z;
qdotdot=BWpara.M\(u-BWpara.C*qdot-RF);
%% Simulation (Newmark)

Idisp=IC(1);Ivel=IC(2);
tol=1e-3;maxit=20;
beta=1/4;gama=1/2; % Newmark-beta
rho=1;LL=1;T=1;FI=1;
coef=[BWpara.uy,BWpara.K,BWpara.alpha,BWpara.n,BWpara.beta,BWpara.gama...
      ,BWpara.deltaA,BWpara.deltaV,BWpara.deltaNu...
      ,BWpara.pw,BWpara.qw,BWpara.zw,BWpara.Fi,BWpara.deltaFi,BWpara.landa];

disp('Newmark')
tic
[ydd,yd,y,Znm,E]=NewmarkSINA(Time,BWpara.M,BWpara.C,BWpara.K,u,Idisp,Ivel,gama,beta,tol,maxit,coef,@NLexpression,rho,LL,T,FI);
toc
RFnm=BWpara.alpha*BWpara.K*y+(1-BWpara.alpha)*BWpara.K*BWpara.uy*Znm;
%% plots
figure,
subplot(411),plot(Time,u);ylabel('F (N)');set(gca,'FontSize',12),
subplot(412),plot(Time,q,'b',Time,y,'r');ylabel('Disp. (m)');set(gca,'FontSize',12),
subplot(413),plot(Time,qdot,'b',Time,yd,'r');ylabel('Vel. (m/s)');set(gca,'FontSize',12),
subplot(414),plot(Time,qdotdot,'b',Time,ydd,'r');ylabel('Acc. (m/s^2)');set(gca,'FontSize',12),
legend('ODE','Newmark')
figure,
subplot(211),plot(Time,z,'b',Time,Znm,'r');xlabel('Time (sec)');ylabel('Z');set(gca,'FontSize',12),
subplot(212),plot(Time,e,'b',Time,E,'r');xlabel('Time (sec)');ylabel('E (J)');set(gca,'FontSize',12),
legend('ODE','Newmark')
figure,
subplot(121),plot(q,z,'b',y,Znm,'r');xlabel('Disp. (m)');ylabel('Z (N)');set(gca,'FontSize',12),
subplot(122),plot(q,RF,'b',y,RFnm,'r');xlabel('Disp. (m)');ylabel('Restoring Force (N)');set(gca,'FontSize',12),
legend('ODE','Newmark')