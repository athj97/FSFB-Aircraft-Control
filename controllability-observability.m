clear; clc; close all;
syms 's'
A = [-.313 56.7 0; -.139 -.426 0; 0 56.7 0];
B = [.232; .0203; 0];
C = [0 0 1];
D = 0;
sys = ss(A, B, C, D);
x0 = [1; 1; 1];

e = eig(A)
%% Hand Calculations
TIME_DURATION = 10;
t = linspace(0,TIME_DURATION);

eAto0 = x0(1).*exp(-(739.*t)/2000).*(cos((31512431^(1/2).*t)/2000)+(113*31512431^(1/2).*sin((31512431^(1/2).*t)/2000))/31512431) + x0(2).*(113400.*31512431^(1/2).*exp(-(739..*t)/2000).*sin((31512431^(1/2).*t)/2000))/31512431 + x0(3).*0;
eAto1 = x0(1).*-(278.*31512431^(1/2).*exp(-(739.*t)/2000).*sin((31512431^(1/2).*t)/2000))/31512431 + x0(2).*exp(-(739.*t)/2000).*(cos((31512431^(1/2).*t)/2000)-(113.*31512431^(1/2).*sin((31512431^(1/2).*t)/2000))/31512431) + x0(3).*0;
eAto2 = x0(1).*((1313550.*exp(-(739.*t)/2000).*(cos((31512431^(1/2).*t)/2000)+(739.*31512431^(1/2).*sin((31512431^(1/2).*t)/2000))/31512431))/1335773-1313550/1335773) + x0(2).*(2957850/1335773-(2957850.*exp(-(739.*t)/2000).*(cos((31512431^(1/2).*t)/2000)-(15797969.*31512431^(1/2).*sin((31512431^(1/2).*t)/2000))/9863390903))/1335773) + x0(3).*1;

%% Zero Input
[i ii initial_x] = lsim(sys,zeros(100,1),t, x0);
figure; subplot(3,1,1); hold on; grid on;
title('Zero Input Response')
plot(ii,initial_x(:,1),'linewidth',2)
plot(t,eAto0,'r--','linewidth',2); 
legend('MATLAB','Hand Calculation');
xlabel('Time [s]')
ylabel('[x1]')
subplot(3,1,2); hold on; grid on;
title('Zero Input Response')
plot(ii,initial_x(:,2),'linewidth',2)
plot(t,eAto1,'r--','linewidth',2); 
legend('MATLAB','Hand Calculation');
xlabel('Time [s]')
ylabel('[x2]')
subplot(3,1,3); hold on; grid on;
title('Zero Input Response')
plot(ii,initial_x(:,3),'linewidth',2)
plot(t,eAto2,'r--','linewidth',2); 
legend('MATLAB','Hand Calculation');
xlabel('Time [s]')
ylabel('[x3]')

%% Zero State / Step
% Input of system is 1 zero_state and step are equal
clear t
syms 't'
syms 'tau'
eAtauo0 = B(1).*exp(-(739.*(t-tau))/2000).*(cos((31512431^(1/2).*(t-tau))/2000)+(113*31512431^(1/2).*sin((31512431^(1/2).*(t-tau))/2000))/31512431) + B(2).*(113400.*31512431^(1/2).*exp(-(739.*(t-tau))/2000).*sin((31512431^(1/2).*(t-tau))/2000))/31512431 + B(3).*0
eAtauo1 = B(1).*-(278.*31512431^(1/2).*exp(-(739.*(t-tau))/2000).*sin((31512431^(1/2).*(t-tau))/2000))/31512431 + B(2).*exp(-(739.*(t-tau))/2000).*(cos((31512431^(1/2).*(t-tau))/2000)-(113.*31512431^(1/2).*sin((31512431^(1/2).*(t-tau))/2000))/31512431) + B(3).*0;
eAtauo2 = B(1).*((1313550.*exp(-(739.*(t-tau))/2000).*(cos((31512431^(1/2).*(t-tau))/2000)+(739.*31512431^(1/2).*sin((31512431^(1/2).*(t-tau))/2000))/31512431))/1335773-1313550/1335773) + B(2).*(2957850/1335773-(2957850.*exp(-(739.*(t-tau))/2000).*(cos((31512431^(1/2).*(t-tau))/2000)-(15797969.*31512431^(1/2).*sin((31512431^(1/2).*(t-tau))/2000))/9863390903))/1335773) + B(3).*1;

zs_0 = int(eAtauo0,tau,0,t)
zs_1 = int(eAtauo1,tau,0,t);
zs_2 = int(eAtauo2,tau,0,t);

t = linspace(0,TIME_DURATION);

zs_0_vector = subs(zs_0,t);
zs_1_vector = subs(zs_1,t);
zs_2_vector = subs(zs_2,t);

% Matlab calculations
[z zz zero_state_matlab] = lsim(sys,ones(100,1),t);

figure; subplot(3,1,1); hold on; grid on;
title('Zero State Response')
xlabel('Time [s]')
ylabel('[x1]')
plot(t,zero_state_matlab(:,1),'linewidth',2);
plot(t,zs_0_vector,'r--','linewidth',2);
legend('MATLAB','Hand Calculation');
subplot(3,1,2); hold on; grid on;
title('Zero State Response')
xlabel('Time [s]')
ylabel('[x2]')
plot(t,zero_state_matlab(:,2),'linewidth',2);
plot(t,zs_1_vector,'r--','linewidth',2)
legend('MATLAB','Hand Calculation');
subplot(3,1,3); hold on; grid on;
title('Zero State Response')
xlabel('Time [s]')
ylabel('[x3]')
plot(t,zero_state_matlab(:,3),'linewidth',2);
plot(t,zs_2_vector,'r--','linewidth',2)
legend('MATLAB','Hand Calculation');

%% Complete Response
t = linspace(0,TIME_DURATION);
complete_matlab_1 = initial_x(:,1) + zero_state_matlab(:,1);
complete_hand_1 = eAto0 + zs_0_vector;
complete_matlab_2 = initial_x(:,2) + zero_state_matlab(:,2);
complete_hand_2 = eAto1 + zs_1_vector;
complete_matlab_3 = initial_x(:,3) + zero_state_matlab(:,3);
complete_hand_3 = eAto2 + zs_2_vector;

figure; subplot(3,1,1); hold on; grid on;
plot(t,complete_matlab_1,'linewidth',2);
plot(t,complete_hand_1,'r--','linewidth',2);
title('Complete Response')
xlabel('Time [s]')
ylabel('[x1]')
legend('MATLAB','Hand Calculation')
subplot(3,1,2); hold on; grid on;
plot(t,complete_matlab_2,'linewidth',2);
plot(t,complete_hand_2,'r--','linewidth',2);
title('Complete Response')
xlabel('Time [s]')
ylabel('[x2]')
legend('MATLAB','Hand Calculation')
subplot(3,1,3); hold on; grid on;
plot(t,complete_matlab_3,'linewidth',2);
plot(t,complete_hand_3,'r--','linewidth',2);
title('Complete Response')
xlabel('Time [s]')
ylabel('[x3]')
legend('MATLAB','Hand Calculation')

%% Impulse Response
syms 't'

Impulse_Hand_Calc = ilaplace(C*inv(s*eye(3)-A)*B+D);
t = linspace(0,TIME_DURATION);
impulse_hand_vector = subs(Impulse_Hand_Calc,t);

[q qq qqq] = impulse(sys,t);

figure; hold on; grid on;
title('Impulse Response Hand vs Matlab')
xlabel('Time [s]');
plot(t,qqq(:,3),'linewidth',2)
plot(t,impulse_hand_vector(1:100),'r--','linewidth',2);
legend('MATLAB','Hand Calculations')
