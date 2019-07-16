clear; clc; close all;

%% Intro:
% This script is designed to classify continous waveform in to RN, NRN and
% the mixture of these two.
% User shoud adjust parameters to adopt different data sets.
% Haoran Meng
% Time: 03/06/2019
% Reference: Meng, H., Ben-Zion, Y. and Johnson, C., 2019. Detection of 
%   random noise and anatomy of continuous seismic waveforms in dense 
%   array data near Anza California: Geophysical Journal International, 
%   in revision.

thrd1 = 0.45; thrd2 = 0.15; % user specified thresholds to differentiate RN, NRN & Mix

%% Compute rms of each window to find noise candidates
load('RN_Xcorr.mat');
for i = 1:nn
    rms_wvfm(i,1) = rms(wvfm(i,:));
end

[B, I] = sort(rms_wvfm);

% find noise candidate
RN_I0 = sort(I(1:1000));
xcc_0 = xcc(RN_I0,RN_I0);
xcc_v0 = xcc_0(:);

for i = 1:length(RN_I0)
    temp = xcc_0(i,:); temp(i) = [];
    var_xcc_0(i,1) = var(temp);         % variance
end

j = 0;
threshold = median(median(xcc_v0));     % median
for i = 1:length(RN_I0)                 % remove outliers
    temp = xcc_0(i,:); temp(i) = [];
    if (var_xcc_0(i) <= 1.1*median(var_xcc_0)) && (median(temp) >= 0.9*threshold) && (median(temp) <= 1.1*threshold)
        j = j + 1;
        RN_I1(j,1) = RN_I0(i);
    end
end
xcc_1 = xcc(RN_I1,RN_I1);

A0 = zeros(1,L/2);
for i = 1:length(RN_I1)
    A0 = A0 + abs(wvfm_F(RN_I1(i),:));
end
A0 = A0/i;
figure(1)
plot([0 f],[0 A0],'k','linewidth',2); hold on; % Initial mean RN spectra

%% Classification
thrd10 = 0.6; thrd20 = 0.1;                 % Initial threshold (more strict)
for i = 1:nn
    A = abs(wvfm_F(i,:));
    d_spec_L2(i,1) = norm((A - A0));        % Spectral deviation
    xcc_mdn(i,1) = median(xcc(i,RN_I1));    % Mean XCC 
end
x = d_spec_L2; y = xcc_mdn;
dx = std(x)/5; dy = std(y)/5;
for i = 1:length(x)                         % compute data points density
    xInRange = (x>=x(i)-dx) & (x<=x(i)+ dx);
    yInRange = (y>=y(i)-dy) & (y<=y(i)+dy);
    xyInRange = xInRange & yInRange;
    num(i) = sum(xyInRange);
end
num = num/max(num);

% First iteration
xyInRange = num >= thrd10;
A0 = zeros(1,L/2);
j = 0;
for i = 1:nn
    if xyInRange(i)
        j = j + 1;
        RN_I2(j) = i;
        A0 = A0 + abs(wvfm_F(i,:));
    end
end
A0 = A0/j;
figure(1)
plot([0 f],[0 A0],'b','linewidth',2); hold on;

xyInRange = num <= thrd20;
j = 0;
for i = 1:nn
    if xyInRange(i)
        j = j + 1;
        NRN_I2(j) = i;
    end
end

for i = 1:nn
    A = abs(wvfm_F(i,:));
    d_spec_L2(i,1) = norm((A - A0));
    xcc_mdn(i,1) = median(xcc(i,RN_I2));
    temp = xcc(i,NRN_I2); [M,I] = max(temp); temp(I) = [];
    xcc_std(i,1) = std(temp);
end
x = d_spec_L2; y = xcc_mdn;
dx = std(x)/5; dy = std(y)/5;
for i = 1:length(x)
    xInRange = (x>=x(i)-dx) & (x<=x(i)+ dx);
    yInRange = (y>=y(i)-dy) & (y<=y(i)+dy);
    xyInRange = xInRange & yInRange;
    num(i) = sum(xyInRange);
end
num = num/max(num);
rho = 1./xcc_std; rho = rho'/max(rho); rho = rho.*num; rho = rho'/max(rho);

% Second iteration, using user specified threshold
xyInRange = rho >= thrd1;
A0 = zeros(1,L/2);
j = 0;
for i = 1:nn
    if xyInRange(i)
        j = j + 1;
        RN_I3(j) = i;
        A0 = A0 + abs(wvfm_F(i,:));
    end
end
A0 = A0/j;
figure(1)
plot([0 f],[0 A0],'c','linewidth',2); hold on;

xyInRange = num <= thrd2;
j = 0;
for i = 1:nn
    if xyInRange(i)
        j = j + 1;
        NRN_I3(j) = i;
    end
end

for i = 1:nn
    A = abs(wvfm_F(i,:));
    d_spec_L2(i,1) = norm((A - A0));
    xcc_mdn(i,1) = median(xcc(i,RN_I3));
    temp = xcc(i,NRN_I3); [M,I] = max(temp); temp(I) = [];
    xcc_std(i,1) = std(temp);
end
x = d_spec_L2; y = xcc_mdn;
dx = std(x)/5; dy = std(y)/5;
for i = 1:length(x)
    xInRange = (x>=x(i)-dx) & (x<=x(i)+ dx);
    yInRange = (y>=y(i)-dy) & (y<=y(i)+dy);
    xyInRange = xInRange & yInRange;
    num(i) = sum(xyInRange);
end
num = num/max(num);
rho = 1./xcc_std; rho = rho'/max(rho); rho = rho.*num; rho = rho'/max(rho);

% Third iteration
xyInRange = rho >= thrd1;
A0 = zeros(1,L/2);
j = 0;
for i = 1:nn
    if xyInRange(i)
        j = j + 1;
        RN_I4(j) = i;
        A0 = A0 + abs(wvfm_F(i,:));
    end
end
A0 = A0/j;
figure(1)
plot([0 f],[0 A0],'m','linewidth',2); hold on;

xyInRange = rho <= thrd2;
j = 0;
for i = 1:nn
    if xyInRange(i)
        j = j + 1;
        NRN_I4(j) = i;
    end
end

for i = 1:nn
    A = abs(wvfm_F(i,:));
    d_spec_L2(i,1) = norm((A - A0));
    xcc_mdn(i,1) = median(xcc(i,RN_I4));
    temp = xcc(i,NRN_I4); [M,I] = max(temp); temp(I) = [];
    xcc_std(i,1) = std(temp);
end
x = d_spec_L2; y = xcc_mdn;
dx = std(x)/5; dy = std(y)/5;
for i = 1:length(x)
    xInRange = (x>=x(i)-dx) & (x<=x(i)+ dx);
    yInRange = (y>=y(i)-dy) & (y<=y(i)+dy);
    xyInRange = xInRange & yInRange;
    num(i) = sum(xyInRange);
end
num = num/max(num);
rho = 1./xcc_std; rho = rho'/max(rho); rho = rho.*num; rho = rho'/max(rho);

% Forth iteration
xyInRange = rho >= thrd1;
A0 = zeros(1,L/2);
j = 0;
for i = 1:nn
    if xyInRange(i)
        j = j + 1;
        RN_I5(j) = i;
        A0 = A0 + abs(wvfm_F(i,:));
    end
end
A0 = A0/j;
xyInRange = rho <= thrd2;
j = 0;
for i = 1:nn
    if xyInRange(i)
        j = j + 1;
        NRN_I5(j) = i;
    end
end
xyInRange = (rho > thrd2) & (rho < thrd1);
j = 0;
for i = 1:nn
    if xyInRange(i)
        j = j + 1;
        MIX_I5(j) = i;
    end
end

%% figures
% spectra
figure(1)
plot([0 f],[0 A0],'r','linewidth',2); hold on;
xlabel('f (Hz)')
ylabel('Amplitude')
legend('Initial','1st','2nd','3rd','4th')
set(gca,'FontSize',15)
xlim([0 250])


% Classification
figure(2)
subplot(2,2,1)
scatter(x,y,10,log10(num),'filled'); hold on;
xlim([0 5*std(x)]); ylim([min(y) max(y)])
box on; colorbar; caxis([-2 0]); colormap(jet); axis square;
xlabel('Spectral Deviation')
ylabel('Median Xcorr. Coeff.')
title('Density')
set(gca,'FontSize',15)

subplot(2,2,2)
scatter(x,y,10,xcc_std,'filled'); hold on;
xlim([0 5*std(x)]); ylim([min(y) max(y)])
box on; colorbar; caxis([min(xcc_std) max(xcc_std)/1.5]); colormap(jet); axis square;
xlabel('Spectral Deviation')
ylabel('Median Xcorr. Coeff.')
title('STD')
set(gca,'FontSize',15)

subplot(2,2,3)
scatter(x,y,10,log10(rho),'filled'); hold on;
xlim([0 5*std(x)]); ylim([min(y) max(y)])
box on; colorbar; caxis([-2 0]); colormap(jet); axis square;
xlabel('Spectral Deviation')
ylabel('Median Xcorr. Coeff.')
title('Weighted Density')
set(gca,'FontSize',15)

% results
figure(3)
tt1 = 950;
tt2 = 1150;
t = 1:nn*Fs*win_len; t = t * delta;
pos1 = [0.1 0.74 0.85 0.16];
pos2 = [0.1 0.58 0.85 0.16];
pos3 = [0.1 0.42 0.85 0.16];
pos4 = [0.1 0.26 0.85 0.16];
pos5 = [0.1 0.1 0.85 0.16];

subplot('Position',pos1)
plot(t,tar_wvfm,'b'); hold on;
set(gca,'xticklabel',[])
xlim([0 t(end)])
%xlim([tt1 tt2])

sig_idx = sort(RN_I5);
sig = []; t_sig = [];
for j = 1:length(sig_idx)
    nst1 = (sig_idx(j)-1)*Fs*win_len + 1;
    nst2 = sig_idx(j)*Fs*win_len;
    temp = tar_wvfm(nst1:nst2);
    sig = horzcat(sig,temp');
    sig(end) = nan;
    t = nst1:nst2; t = t*delta;
    t_sig = horzcat(t_sig,t);
end
plot(t_sig,sig,'r'); hold on;

sig_idx = sort(NRN_I5);
sig = []; t_sig = [];
for j = 1:length(sig_idx)
    nst1 = (sig_idx(j)-1)*Fs*win_len + 1;
    nst2 = sig_idx(j)*Fs*win_len;
    temp = tar_wvfm(nst1:nst2);
    sig = horzcat(sig,temp');
    sig(end) = nan;
    t = nst1:nst2; t = t*delta;
    t_sig = horzcat(t_sig,t);
end
plot(t_sig,sig,'k'); hold on;
set(gca,'FontSize',12)
ylim([-10*std(tar_wvfm) 10*std(tar_wvfm)]);
ylabel('Amplitude (count)','fontsize',12)

pctg_RN = length(RN_I5)/nn*100;
pctg_NRN = length(NRN_I5)/nn*100;
pctg_mix = 100 - pctg_RN -pctg_NRN;

subplot('Position',pos2)
[S,F,T] = spectrogram(tar_wvfm,1024,512,1024,1/delta,'yaxis');
spec = 10*log10(abs(S));
pcolor(T,F,spec); shading interp; hold on;
colormap(jet); caxis([30 45]);
set(gca,'xticklabel',[])
set(gca,'FontSize',12)
box on; ylabel('F (Hz)','fontsize',12); xlim([0 3600]); 
xlabel('t (s)','fontsize',12)
xlim([0 t(end)])

subplot('Position',pos1)
title({['Random Noise(red): ' num2str(pctg_RN,'%3.1f') ...
    '%; Not RN(black): ' num2str(pctg_NRN,'%3.1f')...
    '%; Mix(blue): ' num2str(pctg_mix,'%3.1f') '%']},'fontsize',18)

subplot('Position',pos3)
t = 1:3600;
plot(t,y,'k');
set(gca,'xticklabel',[])
xlim([0 3600]); ylim([min(y) max(y)]);
ylabel('Median Xcor. Coef.','fontsize',12)
set(gca,'FontSize',12)
xlim([0 t(end)])

subplot('Position',pos4)
t = 1:3600;
semilogy(t,x,'k'); ylim([min(x) max(x)]);
set(gca,'xticklabel',[])
ylabel('Spectral Deviation','fontsize',12)
set(gca,'FontSize',12)
xlim([0 t(end)])

subplot('Position',pos5)
t = 1:3600;
plot(t,xcc_std,'k'); ylim([min(xcc_std) max(xcc_std)]);
set(gca,'FontSize',12)
ylabel('STD','fontsize',12)
xlabel('t (s)','fontsize',15)
xlim([0 t(end)])


% histogram
figure(2);
xbin = 0:0.02:1;
subplot(2,2,4)

xcc3 = xcc(NRN_I5,NRN_I5); xv3 = xcc3(:);
[n, xout] = hist(xv3,xbin);
bar(xout, n, 'barwidth', 1, 'edgecolor', 'k','facecolor',[0.5 0.5 0.5],'linewidth',1); hold on;
set(gca,'FontSize',12)


xcc3 = xcc(MIX_I5,MIX_I5); xv3 = xcc3(:);
[n, xout] = hist(xv3,xbin);
bar(xout, n, 'barwidth', 1, 'edgecolor', 'b','facecolor',[0.5 0.5 1],'linewidth',1); hold on;
axis square;
set(gca,'FontSize',12)

xcc3 = xcc(RN_I5,RN_I5); xv3 = xcc3(:);
[n, xout] = hist(xv3,xbin);
bar(xout, n, 'barwidth', 1, 'edgecolor', 'r','facecolor',[1 0.5 0.5],'linewidth',1); hold on;
axis square;
set(gca, 'YScale', 'log')
xlabel('Xcorr. Coeff.')
ylabel('Number')
xlim([0 1])
colorbar; colormap(jet);
legend('Not RN','Mix','Random Noise')

set(gca,'FontSize',15)


















