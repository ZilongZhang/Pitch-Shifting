function [Iden,x_out] = RT_PitchShift(x_in,ratio,step)
%通过scaling原有基频合成变掉后的语音
% x : 分割后的单帧声音信号 double n x 1 的列向量 
%fs : 采样频率 int 1x1
%ratio:缩放比例，ratio>1则声音变得尖锐（音调变高），<1则声音变得低沉（音调变低） double 1x1
%step:合成时的帧移，一般使用帧长的1/4. double 1x1

% step           = 512;    % synthesis step [samples]
% ratio    = 0.8    % pitch-shifting ratio
% s_win        = 2048;   % analysis window length [samples]
% [x,FS] = audioread('x1.wav');
% x = x(1:2048*2);
%----- initialize windows, arrays, etc -----
% s_win = length(x_in);
s_win = 2048;
n1       = round(step / ratio);      % analysis step [samples]
tstretch_ratio = step/n1;
w1       = hanning(s_win, 'periodic'); % analysis window
w2       = w1;    % synthesis window
L        = length(x_in);
x_in  = [zeros(s_win, 1); x_in;zeros(s_win-mod(L,n1),1)] / max(abs(x_in));
x_out = zeros(length(x_in),1);
omega    = 2*pi*n1*[0:s_win-1]'/s_win;
phi0     = zeros(s_win,1);
psi      = zeros(s_win,1);
Iden = zeros(length(x_in),1);

%----- for linear interpolation of a grain of length s_win -----
lx   = floor(s_win*n1/step);
x    = 1 + (0:lx-1)'*s_win/lx;
ix   = floor(x);
ix1  = ix + 1;
dx   = x - ix;
dx1  = 1 - dx;

%降噪
coef = 0.1;
hs_win   = 1;


%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU

 pin  = 0;
pout = 0;
pend = length(x_in)-max(s_win,lx);

while pin<pend
  grain = x_in(pin+1:pin+s_win).* w1;
%   grain = x_in(pin+1:pin+s_win);
%===========================================
  f     = fft(fftshift(grain));
  r     = abs(f)/hs_win;
  f =  f.*r./ (r+coef);
  phi   = angle(f);
  %---- computing phase increment ----
  delta_phi = omega + princarg(phi-phi0-omega);
  phi0  = phi;
  psi   = princarg(psi+delta_phi*tstretch_ratio);
  %---- synthesizing time scaled grain ----
  ft    = (r.* exp(i*psi));
  grain = fftshift(real(ifft(ft))).*w2; %重要改动！！！！！！！！！！！去掉了窗
% grain = fftshift(real(ifft(ft)));
  %----- interpolating grain -----
  grain2 = [grain;0];
  grain3 = grain2(ix).*dx1+grain2(ix1).*dx;
%  plot(grain);drawnow;
% ===========================================
  Iden(pout+1:pout+lx)= Iden(pout+1:pout+lx) + hanning(length(grain3)); 
  x_out(pout+1:pout+lx) = x_out(pout+1:pout+lx) + grain3;
  pin    = pin + n1;
  pout   = pout + n1;
  end
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU


%----- listening and saving the output -----
% x  = x(s_win+1:s_win+L);
Iden = Iden(s_win+1:s_win+L); 
x_out = x_out(s_win+1:s_win+L);
% soundsc(x_out, fs);
% t = (length(x_out)-500) /fs;
% pause(t);
