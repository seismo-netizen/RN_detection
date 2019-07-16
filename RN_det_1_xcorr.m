clear; clc; close all;

%% Intro:
% This script is designed to cut the continous waveform in to small
% segments and compute the xcorr of each window pair
% Haoran Meng
% Time: 03/06/2019


%% data preparation
% Target wvfm must be filtered by "highpass corner 2/win_len n 2 p 1"
% User should modify this script to read in waveform
fnm = 'wvfm.sac';       % waveform file name
hr = 1;                 % total data length in hours
win_len = 1;            % window length in second
nn = hr*3600/win_len;

[hdr, data] = load_sac(fnm);    % load waveform
delta = hdr.delta;              % smapling period
Fs = round(1/delta);            % sampling rate
L = win_len * Fs;               % samples per window
tar_wvfm = data(1:hr*3600*Fs);  % target waveform
clearvars hdr data;

%% X-correlation
tic;

tpr = tukeywin(L,0.025); tpr = tpr';    % taper window
f = Fs*(1:(L/2))/L;                     % frequency

% Compute the amplitude spectra of all windows
for i = 1:nn
    nt1 = ((i - 1)*win_len) * Fs + 1;
    nt2 = (i*win_len) * Fs;
    wvfm(i,:) = tar_wvfm(nt1:nt2);
    
    % fft complex spectra
    Y = fft(wvfm(i,:));
    wvfm_F(i,:) = Y(2:L/2+1)/L*2;
    wvfm_F(i,end) = wvfm_F(i,end)/2;
end

% Xcor. coef.
for i = 1:1:nn
    wnd1 = wvfm(i,:);
    for j = i+1:1:nn
        wnd2 = wvfm(j,:);
        [cr,~] = xcorr(wnd1,wnd2,'coeff',win_len*Fs/2);
        xcc(i,j) = max(abs(cr));
        xcc(j,i) = xcc(i,j);
    end
end

for i = 1:1:nn
    xcc(i,i) = 1;
end

save RN_Xcorr.mat

toc;



















