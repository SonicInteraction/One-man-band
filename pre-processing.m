%% Preprocessing Of Signals
% Filters
first = bp00501(first);
% second = extremeLowPass(second);
%% Compression
close all;
x = second;
thr = 0.2;
gr = 0.8;
for i = 1:length(x);
    if abs(x(i)) > thr;
        x(i) = x(i)*gr;
    end
end
x = x./max(x);
figure(1);
subplot(2,1,1);
plot(x);
subplot(2,1,2)
plot(second)
%% Assign compressed signal
second = x;
