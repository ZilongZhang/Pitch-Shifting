% VX_pitch_pv.m   [DAFXbook, 2nd ed., chapter 7]
%===== This program performs pitch shifting 
%===== using the FFT/IFFT approach
%
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------

clear; clf

%----- user data -----
n2           = 512;    % synthesis step [samples]
pit_ratio    = 0.6    % pitch-shifting ratio
s_win        = 2048;   % analysis window length [samples]
[DAFx_in,FS] = audioread('x1.wav');

% [f1,f2] = CepsFormant(DAFx_in,FS)
% DAFx_in = DAFx_in(2048*20:2048*21);
%----- initialize windows, arrays, etc -----
n1       = round(n2 / pit_ratio);      % analysis step [samples]
tstretch_ratio = n2/n1;
w1       = hanning(s_win, 'periodic'); % analysis window
w2       = w1;    % synthesis window
L        = length(DAFx_in);
DAFx_in  = [zeros(s_win, 1); DAFx_in; ...
   zeros(s_win-mod(L,n1),1)] / max(abs(DAFx_in));
DAFx_out = zeros(length(DAFx_in),1);
iden_total= zeros(length(DAFx_in),1);
omega    = 2*pi*n1*[0:s_win-1]'/s_win;
phi0     = zeros(s_win,1);
psi      = zeros(s_win,1);

%----- for linear interpolation of a grain of length s_win -----
lx   = floor(s_win*n1/n2);
x    = 1 + (0:lx-1)'*s_win/lx;
ix   = floor(x);
ix1  = ix + 1;
dx   = x - ix;
dx1  = 1 - dx;


tic
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
pin  = 0;
pout = 0;
k=0;
pend = length(DAFx_in)-max(s_win,lx);

while pin<pend
  grain = DAFx_in(pin+1:pin+s_win).* w1;
%===========================================
  f     = fft(fftshift(grain));
  r     = abs(f);
  phi   = angle(f);
  %---- computing phase increment ----
  delta_phi = omega + princarg(phi-phi0-omega);
  phi0  = phi;
  psi   = princarg(psi+delta_phi*tstretch_ratio);
  %---- synthesizing time scaled grain ----
  ft    = (r.* exp(i*psi));
  grain = fftshift(real(ifft(ft))).*w2;
  %----- interpolating grain -----
  grain2 = [grain;0];
  grain3 = grain2(ix).*dx1+grain2(ix1).*dx;
%  plot(grain);drawnow;
% ===========================================
  iden = hanning(length(grain3),'periodic');
  iden_total(pout+1:pout+lx) = iden_total(pout+1:pout+lx) + iden;
  DAFx_out(pout+1:pout+lx) = DAFx_out(pout+1:pout+lx) + grain3;
  k=k+1;
  pin    = pin + n1;
  pout   = pout + n1;
%   if k >4
%       sound(DAFx_out((k-4)*n1 + 1:(k-1)*n1),FS);
%       pause(0.0001);
%   end
end
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
toc
% pause();
%----- listening and saving the output -----
DAFx_in  = DAFx_in(s_win+1:s_win+L);
DAFx_out = DAFx_out(s_win+1:s_win+L) / max(abs(DAFx_out));

% [f1,f2] = CepsFormant(DAFx_out,FS)
soundsc(DAFx_out, FS);
% audiowrite('e1 ·ÇÊµÊ±-1.2.wav',DAFx_out,FS);
% wavwrite(DAFx_out, FS, 'flute2_pitch_pv.wav');