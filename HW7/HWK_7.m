clear; close all; clc
data = load('homework7.mat');

t = data.time; y = data.y;

%% Problem 1

figure(1)
periodogram(y)
hold on
title('Raw Series PSD')
dt = t(2) - t(1);

df = 1/dt;

dom_normfreq = [0.636719 0.03125 0.160156];
dom_freq = df/2*dom_normfreq;

disp('The dominant frequencies from most to least:')
disp(dom_freq)

%% Problem 2

dy_dt = diff(y)/dt;
figure(2)
periodogram(dy_dt)
hold on
title('differentiation filter PSD')

% behaves as a high-pass filter, power contribution of lower frequencies is
% reduced while larger frequency power contribution is amplified.

%% Problem 3

% Casual Indexing Running Mean

b = ones(1,9)/9;

y_s = conv(y,b, 'same');

figure(3)
plot(t,y_s)
hold on
title('Running Mean filtered Series')

figure(4)
periodogram(y_s)
hold on 
title('Running Mean filtered PSD')
% Behaves as a low-pass filter

%% Problem 4

% fir1 20th order low pass, cutoff f = 2nd highest f from 1

f_l = fir1(20,2*0.160156/df,"low");

figure(5)
freqz(f_l,1)
hold on
title('Low Pass Filter')


y_l = filter(f_l,1,y);
figure(6)
plot(t,y_l)
hold on
title('Low Pass Filtered Series')


%% Problem 5

% fir1 20th order band pass (3.1804,0.8)
band = 2*[0.16,3.14]./df;
f_b = fir1(20,band);

figure(7)
freqz(f_b,1)
hold on
title('Band Pass Filter')
hold off

figure(8)
y_b = filter(f_b,1,y);
plot(t,y_b)
hold on
title('Band Pass filtered Series')

%% Problem 6

% butter func 20th order low pass cut off 0.1561

[b_l,a] = butter(20, 0.160156*2/df);

figure(9)
freqz(b_l,a,[],df)
hold on 
title('ButterWorth Low Pass')
hold off

figure(10)
y_bl = filter(b_l,a,y);

plot(t,y_bl)
hold on
title('Low Pass Butterworth Filtered Series')
hold off