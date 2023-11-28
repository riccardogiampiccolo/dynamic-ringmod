function [b, v, r] = diodeNRsolver(a, v_guess, Z)

%% Copyright (C) 2020 Riccardo Giampiccolo
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
% Implementation of the Newthon-Raphson solver used in:
%
% A. Bernardini, P. Maffezzoni, L. Daniel and A. Sarti, "Wave-Based Analysis 
% of Large Nonlinear Photovoltaic Arrays," in IEEE Transactions on Circuits 
% and Systems I: Regular Papers, vol. 65, no. 4, pp. 1363-1376, April 2018, 
% doi: 10.1109/TCSI.2017.2756917.

%% Initialization
    eNR = 1e-10;
    
    % Diode Parameters
    Vth = 26e-3; % thermal voltage
    Is = 1e-9; % saturation current
    eta = 2.19; % ideality factor
    Rs = 1e-3; % diode series resistance
    Rp= 100e3; % diode parallel resistance
 
    h = @(v) Is * (exp(((v * (Z + Rs) - a * Rs))/(eta * Vth * Z)) - 1) + (v * (Z + Rp + Rs) - a * (Rp + Rs))/(Z * Rp);
    dh = @(v) (Is * (Z + Rs)/(eta * Vth * Z)) * exp(((v * (Z+ Rs) - a * Rs))/(eta * Vth * Z)) + (Z + Rp + Rs)/(Z * Rp);
    df_i = @(v, i) (1 + Rs/Rp + Is*Rs/(eta*Vth)*exp((v - Rs * i)/(eta * Vth)))/(1/Rp + Is/(eta*Vth) * exp((v - Rs * i)/(eta*Vth)));

%% Algorithm
    
    diff = 1000;
    while diff > eNR
        v = v_guess - h(v_guess)/dh(v_guess); % NR solver
        diff = abs(v - v_guess);
        v_guess = v;
    end
    
%     i = (a - v)/Z;
    b = 2*v - a;
    r = df_i(v, in_new);

end