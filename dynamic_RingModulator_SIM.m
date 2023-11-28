close all; clear; clc;

%% Copyright (C) 2021 Riccardo Giampiccolo
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%
% Implementation of the method presented in:
%
% A. Bernardini, P. Maffezzoni, L. Daniel and A. Sarti, "Wave-Based Analysis 
% of Large Nonlinear Photovoltaic Arrays," in IEEE Transactions on Circuits 
% and Systems I: Regular Papers, vol. 65, no. 4, pp. 1363-1376, April 2018, 
% doi: 10.1109/TCSI.2017.2756917.
%
% and applied to solve the Dynamic Ring Modulator circuit in:
%
% A. Bernardini, P. Maffezzoni and A. Sarti, "Linear Multistep Discretization
% Methods With Variable Step-Size in Nonlinear Wave Digital Structures for
% Virtual Analog Modeling," in IEEE/ACM Transactions on Audio, Speech,
% and Language Processing, vol. 27, no. 11, pp. 1763-1776, Nov. 2019,
% doi: 10.1109/TASLP.2019.2931759.

addpath("LTspice/")
addpath("lib/")

eSIM = 1e-10; % convergence threshold of the iterative method
n_ports = 13;
nl_elements = 4;

fin = 150; % input frequency
fc = 50; % carrier frequency
fs = 44100; % sampling frequency
Ts = 1/fs;
StopTime = 50e-3;
t = 0:Ts:StopTime;

%% Circuit Parameters

Vin = sin(2*pi*fin*t);
Vc = sin(2*pi*fc*t);
Rin = 80;
Rc = 1;
La = 0.8;
Z_La = 2*La*fs; 
Lb = 0.8;
Z_Lb = 2*Lb*fs;
Ca = 1e-9;
Z_Ca = Ts/(2*Ca);
Cb = 1e-9;
Z_Cb = Ts/(2*Cb);
Cd = 1e-9;
Z_Cd = Ts/(2*Cd);
Rd = 50;
Rout = 600; 
Z_D1 = 1; % initialization of diode port resistance
Z_D2 = 1;
Z_D3 = 1;
Z_D4 = 1;
Z = diag([Z_D1, Z_D2, Z_D3, Z_D4, Z_La, Z_Lb, Z_Ca, Z_Cb, Z_Cd, Rin, Rd, Rout, Rc]); % matrix containing all the port resistances

%% Initialization of Vectors

Vout = zeros(length(t),1);
a = zeros(n_ports, 1);
b = zeros(n_ports, 1);

v_diode_guess = 0.2*ones(1, nl_elements);
slope = zeros(1, nl_elements);
v_old_iter = 0;
a_prevSample = zeros(5, 1);

%% Fundamental Cut-Set Matrix

I = eye(4);
F = [0.5, 0.5, -0.5, -0.5, 1, 0, 1, 0, 0; 
       1,  -1,    1,   -1, 0, 0, 0, 0, 0;
    -0.5, 0.5,  0.5, -0.5, 0, 1, 0, 1, 0; 
       1,  -1,    1,   -1, 0, 0, 0, 0, 1];

Q = [F, I]; % fundamental cut-set matrix
% B = [eye(9), -F.']; % fundamental loop matrix

% if numbers of turns are provided:
%
% gamma = 4; % example
% theta = 4;
% tau = 8;
% xi = 8;
% lambda = 16;
% mu = 8;
% 
% B = [mu*lambda 0 0 0 0 0 0 0 0 -theta*lambda -mu*lambda xi*mu -mu*lambda; 
%      0 mu*lambda 0 0 0 0 0 0 0 -gamma*lambda  mu*lambda -xi*mu  mu*lambda;
%      0 0 mu*lambda 0 0 0 0 0 0  gamma*lambda -mu*lambda -tau*mu -mu*lambda;
%      0 0 0 mu*lambda 0 0 0 0 0  theta*lambda  mu*lambda  tau*mu  mu*lambda;
%      0 0 0 0 mu*lambda 0 0 0 0  -mu*lambda   0   0   0;
%      0 0 0 0 0 mu*lambda 0 0 0   0   0  -mu*lambda   0;
%      0 0 0 0 0 0 mu*lambda 0 0  -mu*lambda   0   0   0;
%      0 0 0 0 0 0 0 mu*lambda 0   0   0  -mu*lambda   0;
%      0 0 0 0 0 0 0 0 mu*lambda   0  -mu*lambda   0   0];
% I_nt = mu*lambda*eye(4);
% Q = [-B(:, 10:end).', I_nt];


%% Scattering Matrix

S = @(Z) 2*Q'*((Q/Z*Q')\(Q/Z)) - eye(n_ports); 
% S = @(Z) eye(n_ports) - 2*Z*B'*((B*Z*B')\B);

%% Algorithm

for n = 1 : length(t)
    
    S_temp = S(Z);  

    % (Linear) Local Scattering Stage
    
    b(5) = -a_prevSample(1);
    b(6) = -a_prevSample(2);
    
    b(7) = a_prevSample(3);
    b(8) = a_prevSample(4);
    b(9) = a_prevSample(5);

    b(10) = Vin(n);
    b(13) = Vc(n);
    
    flag = 1;
    while flag

        % (Nonlinear) Local Scattering Stage
        
        for x = 1 : nl_elements
%             [b(x), v_diode_guess(x), slope(x)] = diodeNRsolver(a(x), v_diode_guess(x), Z(x,x)); % Newton-Raphson solver
            [b(x), slope(x)] = DiodeWrightOmega(a(x), Z(x,x)); % Wright Omega function
        end
        
        % Global Scattering Stage 

        a = S_temp*b; 

        % Convergence Check
        
        v = (a + b)/2;

        if max(abs(v - v_old_iter)) < eSIM
            % currents must be computed before the port resistance updates
            
            for x = 1 : nl_elements
                Z(x, x) = slope(x); % update of port resistances
            end
            
            flag = 0;
        end
        
        v_old_iter = v;
    end
    
    a_prevSample = a(5:9); % update of waves for dynamic elements
    
    Vout(n) = v(12);

end

%% Load DATA from LTspice

loadedFile = load('ltspice_Vout.txt');
SPICE_time = loadedFile(:,1);
SPICE_amp = loadedFile(:,2);

%% Plot

figure
set(gcf,'color','w');
plot(t, Vout, 'b','LineWidth',2,'DisplayName','WD');
hold on
plot(SPICE_time, SPICE_amp, 'r--','LineWidth',3,'DisplayName','LTspice');
xlabel('Time [s]','interpreter','latex','FontSize',18);
ylabel('Voltage [V]','interpreter','latex','FontSize',18);
legend('show','interpreter','latex','FontSize',13);
