clc;
clear;
close all;

[y_talking, fs_talking] = audioread('Talking.m4a');
[y_noise, fs_noise] = audioread('Noise.m4a');
[y_singing, fs_singing] = audioread('Singing.m4a');
[y_male, fs_male] = audioread('406410007.wav');
[y_female, fs_female] = audioread('408410036.wav');
[y_936, fs_936] = audioread('936.m4a');
[y_936_1, fs_936_1] = audioread('936.mp4');



%% time domain

% talking
t_talking = linspace(0, length(y_talking)/fs_talking, length(y_talking));
subplot(311), plot(t_talking, y_talking);
title('Talking Time Domain')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')

% noise
t_noise = linspace(0, length(y_noise)/fs_noise, length(y_noise));
subplot(312), plot(t_noise, y_noise);
title('Noise Time Domain')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')

% singing
t_singing = linspace(0, length(y_singing)/fs_singing, length(y_singing));
subplot(313), plot(t_singing, y_singing);
title('Singing Time Domain')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')

%% spectral view

% x, window, noverlap, nfft, fs
figure;
subplot(311), spectrogram(y_talking, 1024, [], [], fs_talking, 'yaxis');
title('Talking Spectrogram');
subplot(312), spectrogram(y_noise, 1024, [], [], fs_noise, 'yaxis');
title('Noise Spectrogram');
subplot(313), spectrogram(y_singing, 1024, [], [], fs_singing, 'yaxis');
title('Singing Spectrogram');

%% freq domain

figure;
% talking
f_talking = linspace(0, 660, 1024);
Y_talking = abs(fft(y_talking, 1024));
subplot(311), plot(f_talking(1:512), fftshift(Y_talking(1:512)));
title('Talking Freq Domain');
xlabel('Freq (Hz)');
ylabel('Amplitude (a.u)');

% Noise
f_noise = linspace(0, 660, 1024);
Y_noise = abs(fft(y_noise, 1024));
subplot(312), plot(f_noise(1:512), fftshift(Y_noise(1:512)));
title('Noise Freq Domain');
xlabel('Freq (Hz)');
ylabel('Amplitude (a.u)');

%Singing
f_singing = linspace(0, 660, 1024);
Y_singing = abs(fft(y_singing, 1024));
subplot(313), plot(f_singing(1:512), fftshift(Y_singing(1:512)));
title('Singing Freq Domain');
xlabel('Freq (Hz)');
ylabel('Amplitude (a.u)');

Y_male = fftshift(abs(fft(y_male, 1024)));
Y_female = fftshift(abs(fft(y_female, 1024)));


%% EMD

me_male = ws_distance(fftshift(Y_singing), Y_male)
me_female = ws_distance(fftshift(Y_singing), Y_female)

%% 936-clean your mind
figure;
f_936_1 = linspace(0, fs_936_1, 1024);
Y_936_1 = abs(fft(y_936_1, 1024));
subplot(211), plot(f_936_1(1:512), Y_936_1(1:512));
title('936 Freq Domain');
xlabel('Freq (Hz)');
ylabel('Amplitude (a.u)');

subplot(212), spectrogram(y_936, 1024, [], [], fs_936, 'yaxis');
title('936 Spectrogram');

%%
function wsd = ws_distance(u_samples, v_samples, p)
% WS_DISTANCE 1- and 2- Wasserstein distance between two discrete 
% probability measures 
%   
%   wsd = WS_DISTANCE(u_samples, v_samples) returns the 1-Wasserstein 
%   distance between the discrete probability measures u and v 
%   corresponding to the sample vectors u_samples and v_samples
%
%   wsd = WS_DISTANCE(u_samples, v_samples, p) returns the p-Wasserstein 
%   distance between the discrete probability measures u and v
%   corresponding to the sample vectors u_samples and v_samples. 
%   p must be 1 or 2.
%
% from https://github.com/nklb/wasserstein-distance

if ~exist('p', 'var')
    p = 1;
end

u_samples_sorted = sort(u_samples(:));
v_samples_sorted = sort(v_samples(:));

if p == 1
    
    all_samples = unique([u_samples_sorted; v_samples_sorted], 'sorted');
    
    u_cdf = find_interval(u_samples_sorted, all_samples(1:end-1)) ...
        / numel(u_samples);
    v_cdf = find_interval(v_samples_sorted, all_samples(1:end-1)) ...
        / numel(v_samples);
    
    wsd = sum(abs(u_cdf - v_cdf) .* diff(all_samples));
    
elseif p == 2
    
    u_N = numel(u_samples);
    v_N = numel(v_samples);    
    all_prob = unique([(0:u_N) / u_N, (0:v_N) / v_N], 'sorted').';
    
    u_icdf = u_samples_sorted(fix(all_prob(1:end-1) * u_N) + 1);
    v_icdf = v_samples_sorted(fix(all_prob(1:end-1) * v_N) + 1);
    
    wsd = sqrt(sum((u_icdf-v_icdf).^2 .* diff(all_prob)));
    
else
    
    error('Only p=1 or p=2 allowed.')
    
end
end

function idx = find_interval(bounds, vals)
% Given the two sorted arrays bounds and vals, the function 
% idx = FIND_INTERVAL(bounds, vals) identifies for each vals(i) the index 
% idx(i) s.t. bounds(idx(i)) <= vals(i) < bounds(idx(i) + 1).

m = 0;
bounds = [bounds(:); inf];
idx = zeros(numel(vals), 1);

for i = 1:numel(vals)
    while bounds(m+1) <= vals(i)
        m = m + 1;
    end
    idx(i) = m;
end
end
