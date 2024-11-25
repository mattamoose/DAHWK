clear; clc; close all;

bulkdat = load('slevel.mat');
t = bulkdat.time; % year decimal form
sl = bulkdat.sealevel; % mm

mei = load('mei.mat');
tm = mei.time; % year
m = mei.meival;
%% Problem 1
%pmtm, 3.5 for time-bandwidth
figure(1)
a = polyfit(t,sl,1);
plot(t,sl)
hold on
% a(1)*t.^2
% plot(t,(a(1)*t.^2+a(1)*t + a(2)))
hold off
sld = sl - (a(1)*t + a(2)); % detrended sea level
plot(t,sld)
p = chi2gof(sld);
if p == 0
    disp('From chi2gof, the error in the trend is normally distributed.')
end

figure(2)
tiledlayout (2,1)
nexttile
[px,w1] = pmtm(sld,3.5,'ConfidenceLevel',0.95);
plot(w1,px)
pmtm(sld,3.5,'ConfidenceLevel',0.95);
%xscale('log')
nexttile
[pxx,w2] = periodogram(sld);
plot(w2,pxx)
%xscale('log')
periodogram(sld);

%% Problem 2

figure(3)
spectrogram(sld)
title('Spectrogram of Sea Level')
hold on
xline(0.057)

%% Problem 3
figure(4)
cwt(sld)
yline(0.04)
%% Problem 4

figure(5)
a = polyfit(tm,m,1);
plot(tm,m)

trend = (a(1)*tm + a(2));
md = m - (a(1)*tm + a(2))';
v = chi2gof(sld);
if v == 0
    disp('From chi2gof, the error in the trend is normally distributed.')
end

figure(6)
tiledlayout (2,1)
nexttile
pmtm(md,3.5)

nexttile
periodogram(md)

%% Problem 5
figure(7)

spectrogram(md)
hold on
xline(0.04)
title('Spectrogram of MEI')
figure(8)
cwt(md)
hold on
yline(0.04)

%% Problem 7

A = [ones(size(t)) t cos(2*pi*t) sin(2*pi*t)];

y = A\sl;

sld = sl - A*y;
m = m(524:end);
tm = tm(524:end);

meif = interp1(tm,m,t);

sp = (meif(515) - meif(514))/(t(515) - t(514));
bp = meif(514) - sp*t(514);

meif(516:538) = sp*t(516:538) + bp;

figure(9)
cpsd(meif,sld)

%mod