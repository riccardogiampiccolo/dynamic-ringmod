function [b, r] = DiodeWrightOmega(a, Z)

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
    % Implementation of the Wright Omega function presented in:
    %
    % R. Giampiccolo, M. G. d. Bari, A. Bernardini and A. Sarti, "Wave Digital
    % Modeling and Implementation of Nonlinear Audio Circuits With Nullors," 
    % in IEEE/ACM Transactions on Audio, Speech, and Language Processing,
    % vol. 29, pp. 3267-3279, 2021, doi: 10.1109/TASLP.2021.3120627.
    
    Vt = 26e-3; % thermal voltage
    Is = 1e-9; % saturation current
    eta = 2.19; % ideality factor
    Rs = 1e-3; % series resistance
    Rp = 1e5; % parallel resistance
    
    alpha = 1 - Rs/Z;
    beta = 1 + Rs/Z;
    delta = 1/Z;

    aa = beta/(2*eta*Vt);
    bb = a*alpha/(2*eta*Vt);
    cc = (-0.5*(beta/Rp) - 0.5*(delta))/Is;
    dd = (a*(0.5*delta - 0.5*(alpha/Rp))/Is) + 1;

    arg = bb - aa*(dd/cc) + log(-(aa/cc));
    expArg = exp(arg);

    if expArg==inf
        b = -(dd/cc) - LightweightWrightOmega(arg,4)/aa;
    else
        b = -(dd/cc) - LambertW(expArg)/aa;
    end
    
    
    v = (a + b)/2;
    i = (a - b)/(2*Z);
    
    r = (1 + Rs/Rp + (Is*Rs/(eta*Vt))*exp((v-Rs*i)/(eta*Vt)))/(1/Rp + (Is/(eta*Vt))*exp((v-Rs*i)/(eta*Vt))); % slope
   
end


