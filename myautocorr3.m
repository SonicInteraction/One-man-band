function [h lags blocks energy corrplot] = myautocorr2(x, wlen, fs, plotit,signalname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
h = wlen/2;
win = sqrt(hanning(wlen));
pin = h;
pout = pin + h;
lags = [];
energy = [];
blocks = 0;
if plotit(1) == 1;
    filename = 'autocorr.gif';
end
while pout < length(x);
    grain = x(pin-h+1:pout).*win;
    R = xcorr(grain);
    R2 = R(end/2:end);
    [i v] = findpeaks(R2);
    e = sum(real(fft(grain)).^2);
    energy = [energy var(grain)];
    if isempty(v);
        lags = [lags fs];
    else
        lags = [lags i(1)];
    end
    pin = pin + h;
    pout = pin + h;
    blocks = blocks + 1;
    if plotit(1) == 1;
    figure(1);
    plot(R2);
    hold on;
    if pout < length(x) - wlen;
    text(i+.02,v,num2str((1:numel(v))'))
    plot(i(1),v(1),'X','LineWidth',3);
    end
    title(['Auto Correlation ' signalname]);
    xlabel('lag time in samples');
    ylabel('Correlation');
    hold off;
    pause(.01);
    drawnow;
    % make a gif
     frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if blocks == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.01);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.01);
      end
      
    else
        corrplot = [];
    end
end
energy = energy./max(energy);

end

