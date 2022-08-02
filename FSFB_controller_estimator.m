close all 
clear all
clc

A = [-0.313 56.7 0; -0.0139 -0.426 0; 0 56.7 0];
B = [0.232; 0.0203; 0];
C = [0 0 1];
D = 0;
%% Controller
K = acker(A,B,[-1; -.2; -.3]);
sys1 = ss(A-B*K,B,C,D)
t = [0:0.1:15];

figure(1)
step(sys1,t)
Nbar = -inv(C*inv(A-B*K)*B) %Nbar equation from notes
sys2 = ss(A-B*K,B*Nbar,C,D) %state space with Nbar scaling
hold on
step(sys2,t)
legend('u=r-Kx','u=Nbar*r-Kx','Location','NorthEast')
%% Estimator
x0 = [.1; .1; .2];
xe = [0; 0; 0];
u = ones(151,1);

A = A-B*K;
B = B*Nbar;
L = acker(A',C',[-5 -.5 -1]); L = L';
Aol = [A zeros(size(A)); zeros(size(A)) A];
Bol = [B;B];
Col = [C zeros(size(C)); zeros(size(C)) C];
Dol = zeros(2,1);

Acl = [A zeros(size(A));L*C A-L*C];
Bcl = [B;B];
Ccl = [C zeros(size(C)); zeros(size(C)) C];
Dcl = zeros(2,1);

[yol,xol] = lsim(Aol,Bol,Col,Dol,u,t,[x0;xe]);
[ycl,xcl] = lsim(Acl,Bcl,Ccl,Dcl,u,t,[x0;xe]);

figure(2)
subplot(211)
plot(t,xol(:,[1 2 3]), t,xol(:,[4 5 6]),'--','Linewidth',2);
title('Open loop estimator');ylabel('states');xlabel('Time')
legend('x1','x2','x3','xe1','xe2','xe3')
subplot(212)
plot(t,xol(:,[1 2 3])-xol(:,[4 5 6]),'Linewidth',2)
legend('x1 error','x2 error','x3 error')
ylabel('Estimation Error');xlabel('Time')

figure(3)
subplot(211)
plot(t,xcl(:,[1 2 3]), t,xcl(:,[4 5 6]),'--','Linewidth',2);
title('Closed loop estimator');ylabel('states');xlabel('Time')
legend('x1','x2','x3','xe1','xe2','xe3')
subplot(212)
plot(t,xcl(:,[1 2 3])-xcl(:,[4 5 6]),'Linewidth',2)
legend('x1 error','x2 error','x3 error')
ylabel('Estimation Error');xlabel('Time')
