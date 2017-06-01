clear all
[x,fs] = audioread('e1.wav');
N = 4096;
ratio = 1.2;
step = 512; 
% step = fix(N * 0.75);
num = fix((length(x)-N)/step);
x_pre = 0;
t = 0.03;
% step = fix(N * 0.75);
tic
x_out_total = zeros(1,441000);
Iden_total = zeros(1,441000);
for i = 1:num
    range = (i-1) * step+1:(i-1)*step+N;
    x_in = x(range);
    x_in  = [zeros(step, 1); x_in;zeros(step,1)];
    [Iden,x_out] = RT_PitchShift(x_in,ratio,step);
    x_out = x_out(step+1:end-step);
    x_out = x_out .*hanning(length(x_out), 'periodic')*0.25;
%     x_out = x_out(N+1:end-N)*(step/N);
    x_out = x_out';
%     Iden = Iden(N+1:end-N)*(step/N);
%     Iden = Iden';
% %     x_out = x_out(step+1:end - step);
%     Iden_total(range) = Iden_total(range) + Iden;
    x_out_total(range) = x_out_total(range) + x_out;
%     x_out_total(1:(i+1)*N) = x_out_total(1:(i+1)*N) / max(x_out_total(1:(i+1)*N));
%--------------------------------------------------------------------------
%     sound(x_out_total((i-1)*N+step:(i+1)*N-step),fs);
%     pause(t);
end
sound(x_out_total,fs);
% plot(Iden_total)
toc
audiowrite('e ʵʱ-1.2.wav',x_out_total,fs);