%% Set recording pair
clc;
close all;
clear first firstname second secondname;  fs = fs(1);   close all;
% bass1 bassHarmonic bassPick Occ BassOcc Erhu Flute Gtr Trmp mTrmp gong
% tshaker percussion(1:2.6200*10^5)


% Erhu and bass Harmonic might be persuaded to work
% Flute and bass Harmonic work quite well
% Trmp and mTrmp, my only success
% mTrmp and tshaker, mmmmmmm?
[first firstname second secondname] = WhichFiles(bassHarmonic,'Bass Gtr Pinched Harmonics',Gtr,'Guitar');
fs = fs(1);

% [mm,peak1_ind]=max(pks)
% 'peak value1 at location'
% pks(peak1_ind) %peak
% locs (peak1_ind) %location
%  'peak value2 at location'
% pks(peak1_ind+1)%peask next to the top peak
% locs (peak1_ind+1) %location
% period=locs(peak1_ind+1)-locs(peak1_ind)
% pitch_Hz=fs/period %display pitch in Hz

%% auto correlation
%% Try and normalise r, and add it as a value in the function. Can explain that it seems more robust than straight energy as it also addresses the fidelity of the signal it's tracking.
%% Left to do - modulus, twd for pitch - tidy diagrams
close all;
x = Flute;
[lag I energy] = myautocorr(x,4096,fs);
%%
close all;
l = 1;
threshE = 0.1;
pitch = fs(1)./lag;
fpitch = medfilt1(pitch,l);
%%
close all
l = 15;
pitchTrack = figure(1);
wave = subplot(2,1,1,'Parent',pitchTrack);
hold(wave,'on');
plot(wave,x);
plot(wave,[1:length(pitch);]*2048,fpitch/max(pitch));
legend('Recording','Filtered Pitch Estimation');
pitchplot = subplot(2,1,2,'Parent',pitchTrack);
hold(pitchplot,'on');
plot(pitchplot,pitch);
hold(pitchplot,'on');
plot(pitchplot,max(pitch)*energy/max(energy));
legend('Estimated Pitch','Energy in Signal');
%%


wrappedpitch = medfilt1(mod(pitch2,40),l);
figure(2)
subplot(4,1,1);
plot(Flute); 
xlim([0 length(Flute)]);
title('Sound Source');
subplot(4,1,2);
plot(pitch,'g');

hold on;
plot(medfilt1(pitch,l),'b','LineWidth',2);
legend('Estimated Pitch','Median Filtered');
ylim([0 1000]);
title('Autocorrelated Pitch Estimation');
subplot(4,1,3);
title('Various Thresholds');
plot(energy/max(energy),'r');

plot((I./energy)/max(I./energy));
hold off;
subplot(4,1,4);
%plot(medfilt1(pitch,l),'b','LineWidth',2);
plot(wrappedpitch,'b','LineWidth',2);
%%
[lag2 pitch1] = myautocorr(second,2048,10,fs);

%%
DM = simmx(lag1,lag2);
[m n DX] = dp(DM);
% Create a two colum matrix for time warping
twm = [m' n'];
% strip out identical frames
[dummy i] = unique(twm(:,1));
twm = twm(i,:);
[dummy i] = unique(twm(:,2));
twm = twm(i,:);
