%% Set recording pair
clc;
close all;
clear first firstname second secondname;
% bass1 bassHarmonic bassPick Occ BassOcc Erhu Flute Gtr Trmp mTrmp gong
% tshaker percussion(1:2.6200*10^5) trumpets


% Erhu and bass Harmonic might be persuaded to work
% Flute and bass Harmonic work quite well
% Trmp and mTrmp, my only success
% mTrmp and tshaker, mmmmmmm?
% percussion and trumpets is close
[first, firstname, second, secondname] = WhichFiles(Gtr,'Guitar',Occ,'Occarina');


%% First set the fft values
close all;
clc;

% Parameters
wlen = 1024;   
h = wlen/2;                        % hop size 50%
win = hanning(wlen);


%% Calculate fft for samples
clc;
close all;
[nblocks1, phase1] = fftForTW4(first,win);
[nblocks2, phase2] = fftForTW4(second,win);
%%
% plot phase
clc;
close all;
figure(2);
plot(phase1,'g');
% hold on;
% plot(second,'r');
% title('Amplitude Comparisson Plots');
% xlabel('Time in Samples');
% ylabel('Amplitude');
% legend([firstname],[secondname]);
%% Plot Amplitudes
clc;
close all;
figure(2);
plot(first,'g');
ylim([-0.5 0.5]);
% hold on;
% plot(second,'r');
% % title('Amplitude Comparisson Plots');
% % xlabel('Time in Samples');
% % ylabel('Amplitude');
% legend([firstname],[secondname]);
%% Plot Spectrogram
clc
close all;
figure(20);


% First sample
subplot(2,2,1);
plot(first,'k');
title([firstname]);
xlabel('time in samples');
ylabel('Amplitude');
subplot(2,2,3);
specgram(first);colorbar;
colormap('gray');

set(gca,'xtick',[])
ylabel('Frequency');
% second sample
subplot(2,2,2);
plot(second,'k');
title([secondname]);
xlabel('time in samples');
ylabel('Amplitude');
subplot(2,2,4);
specgram(second);colorbar;
colormap('gray');
set(gca,'xtick',[])
ylabel('Frequency');
%% Calculate the dissimilarity matrix
clc;
close all;
paramVis.fsAudio = fs(1);
[twm m n DX DM] = dissimM(fft1,fft2,(1:100));

%% Plot the Dissimilarity matrix - try and make a gif rotation of this
clc;
close all;
figure(5);
filename = 'rotation.gif';
for n = 1:90
mesh(DM);
colormap('summer');
title(['Disimilarity Matrix of ' firstname ' and ' secondname ' - Phase Only']);
ylabel([firstname]);
xlabel([secondname]);
zlabel('Disimilarity');
view(270+n,60-round((2*abs(n-45))/3));
%drawnow; 

      frame = getframe(5);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.01);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.01);
      end
end
%% Plot cumulative dissimilarity
clc;
close all;
figure(6);
imagesc(DX);
colormap('lines');
%saveas(6,['plots/Cumulative Disimilarity Matrix of ' firstname ' and ' secondname '.png']);
hold on;
plot(n,m,'k','lineWidth',3);
title(['Cumulative Disimilarity Path for ' firstname ' and ' secondname]);
ylabel([firstname]);
xlabel([secondname]);
zlabel('Dissimilarity');
%saveas(6,['plots/Cumulative Disimilarity Matrix of ' firstname ' and ' secondname ' with path.png']);
%% translate this to sample numbers
clc;
close all;
flippedtwm = flip(twm);
cf1 = length(first)/nblocks1;
cf2 = length(second)/nblocks2;
ntwm = [cf1.*twm(:,1) cf2.*twm(:,2)];
%% Time Warping Path
close all;
clc;
visualizeAP(fliplr(floor(ntwm)),fs(1));
title(['Time Warping Path to match ' firstname ' and ' secondname]);
% saveas(1,'plots/Trumpets Time Warping Path.png');
%% Let's do the time warp again!!!
close all;
clc;


z1 = hpTSM(second,fliplr(floor(ntwm)));


%% Plot results
close all;
clc;
figure(7);
subplot(2,1,1);
plot(first,'g');
hold on;
plot(second,'r');

title([firstname ' and ' secondname ' Without Time Warping']);
xlabel('Time in Samples');
ylabel('Amplitude');
%legend([firstname],[secondname]);
subplot(2,1,2);
plot(first,'g');
alpha(1);
hold on;
plot(z1,'r');
title([firstname ' and ' secondname ' after Time Warping']);
xlabel('Time in Samples');
ylabel('Amplitude');
% legend([first],[second]);
% saveas(7,'plots/Plots of trumpet time warp.png');
%% Make a Stereo File of the results
clc;
close all;
stereo = [z1 first];
% and for comparison
test = [first second(1:length(first))];

%% Make a mono file of the results
clc;
close all;
second=vertcat(second,zeros(max(size(first))-max(size(second)),1));
mono = [first+z1]./2;
%% If it works save it now!
clc;
close all;
audiowrite('results\Trumpets.wav',mono,fs(1));
%% and load it as a variable
clc;
close all;
trumpets = audioread('results\Trumpets.wav');