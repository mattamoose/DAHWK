clear all; clc

%import data structure and extract dataspace
datstruct_opd= load('orbital_parameter_data.mat');
kyear = datstruct_opd.m(:,1);
oecc = datstruct_opd.m(:,2);
pre_ang = datstruct_opd.m(:,3);
obl_ang = datstruct_opd.m(:,4);

figure(1) 
plot(kyear,oecc)

figure(2)
plot(kyear,pre_ang)

figure(3)
plot(obl_ang)
title('Raw Data')


fs = 1; %kyr
n = numel(kyear);


ffte = fft(oecc);
fftpa = fft(pre_ang);
fftoa = fft(obl_ang);

Pe2 = abs(ffte).^2/n;
Ppa2 = abs(fftpa).^2/n;
Poa2 = abs(fftoa).^2/n;
Pe = Pe2(1:floor(n/2)+1);
Pe = [Pe2(1);2*Pe(2:end)];
%Ppa = [Ppa2(1); 2*Ppa2(2:4096)];
%Poa = [Poa2(1); 2*Poa2(2:4096)];

f = 0:1/n:1/2;

figure(4)
plot(f,Pe)
title('Orbital Eccentricity PSD')
xlabel('Frequency (cycles/kyear)')
set(gca,'XScale','log','YScale','log')

%% Problem 2

% de-mean the data
xe = polyfit(kyear,oecc,1);
xoa = polyfit(kyear,obl_ang,1);
xpa = polyfit(kyear,pre_ang,1);

oe_trend = xe(2).*ones(size(kyear)) + xe(1).*kyear;
oa_trend = xoa(2).*ones(size(kyear)) + xoa(1).*kyear;
pa_trend = xpa(2).*ones(size(kyear)) + xpa(1).*kyear;
oecc = oecc - oe_trend;
pre_ang = pre_ang - pa_trend;
obl_ang = obl_ang - oa_trend;

ffte = fft(oecc);
fftpa = fft(pre_ang);
fftoa = fft(obl_ang);

Pe2 = abs(ffte).^2/n;
Ppa2 = abs(fftpa).^2/n;
Poa2 = abs(fftoa).^2/n;

Pe = Pe2(1:floor(n/2)+1);
Ppa = Ppa2(1:floor(n/2)+1);
Poa = Poa2(1:floor(n/2)+1);

Pe = [Pe(1);2*Pe(2:end)];
Ppa = [Ppa(1); 2*Ppa(2:end)];
Poa = [Poa(1); 2*Poa(2:end)];

figure(5)
plot(f,Poa)
title('Obliquity Angle PSD')
xlabel('Frequency (cycles/kyr)')
set(gca,'XScale','log','YScale','log')

figure(6)
plot(f,Pe)
title('Orbital Eccentricity PSD')
xlabel('Frequency (cycles/kyear)')
set(gca,'YScale','log','XScale','log')

figure(7)
plot(f, Ppa)
title('Precession Angle PSD')
xlabel('Frequency (cycles/kyear)')
set(gca,"YScale","log", "XScale","log")


% Band Pass
xoa = find(Poa == max(Poa));
xoe = find(Pe == max(Pe));
xpa = find(Ppa == max(Ppa));
fmax_oa = f(xoa(1));
fmax_oe = f(xoe(1));
fmax_pa = f(xpa(1));

i_oa = (f >= fmax_oa - 0.01) & (f <= fmax_oa + 0.01);
i_oe = (f >= fmax_oe -0.01) & (f <= fmax_oe + 0.01);
i_pa = (f >= fmax_pa - 0.01) & (f <= fmax_pa + 0.01);




Poa(i_oa) = zeros(size(Poa(i_oa)));
Pe(i_oe) = zeros(size(Pe(i_oe)));
Ppa(i_pa) = zeros(size(Ppa(i_pa)));

% magnitude spectrum
mag_oa = sqrt(Poa*n);
mag_oe = sqrt(Pe*n);
mag_pa = sqrt(Ppa*n);

% phase info
phase_oa = 2*pi*rand(size(mag_oa));
phase_oe = 2*pi*rand(size(mag_oe));
phase_pa = 2*pi*rand(size(mag_pa));

cs_oa = mag_oa.*exp(1i*phase_oa);
cs_oe = mag_oe.*exp(1i*phase_oe);
cs_pa = mag_pa.*exp(1i*phase_pa);

%[sqrt(Poa(1))*n; sqrt(Poa(2:end)/2)*n; sqrt(fliplr(conj(Poa(2:end)))/2)*n]
Yoa = [sqrt(Poa(1)*n); sqrt(Poa(2:end)/2*n); sqrt(fliplr(conj(Poa(2:end)/2*n)))];
Ye = [sqrt(Pe(1)*n); sqrt(Pe(2:end)/2*n); sqrt(fliplr(conj(Pe(2:end)/2*n)))];
Ypa = [sqrt(Ppa(1)*n); sqrt(Ppa(2:end)/2*n); sqrt(fliplr(conj(Ppa(2:end)/2*n)))];
Oa_rcnstrct = ifft(cs_oa,5001,'symmetric');
Pe_rcnstrct = ifft(cs_oe,5001,'symmetric');
Ppa_rcnstrc = ifft(cs_pa,5001,'symmetric'); 

figure(8)
plot(kyear,Oa_rcnstrct)
title('Reconstructed Obliquity Angle Time Series')
xlabel('time (kyr)')

figure(9)
plot(kyear,Pe_rcnstrct)
title('Reconstructed Orbital Eccentricity Time Series')
xlabel('time (kyr)')

figure(10)
plot(kyear,Ppa_rcnstrc)
title('Reconstructed Precession Angle Time Series')
xlabel('time (kyr)')

%% Reconstructed Periodogram

ffte = fft(Pe_rcnstrct);
fftpa = fft(Ppa_rcnstrc);
fftoa = fft(Oa_rcnstrct);

Pe2 = abs(ffte).^2/n;
Ppa2 = abs(fftpa).^2/n;
Poa2 = abs(fftoa).^2/n;

Pe = Pe2(1:floor(n/2)+1);
Ppa = Ppa2(1:floor(n/2)+1);
Poa = Poa2(1:floor(n/2)+1);

Pe = [Pe(1);2*Pe(2:end)];
Ppa = [Ppa(1); 2*Ppa(2:end)];
Poa = [Poa(1); 2*Poa(2:end)];

figure(11)
plot(f,Poa)
title('Obliquity Angle PSD')
xlabel('Frequency (cycles/kyr)')
set(gca,'XScale','log','YScale','log')

figure(12)
plot(f,Pe)
title('Orbital Eccentricity PSD')
xlabel('Frequency (cycles/kyear)')
set(gca,'YScale','log','XScale','log')

figure(13)
plot(f, Ppa)
title('Precession Angle PSD')
xlabel('Frequency (cycles/kyear)')
set(gca,'YScale','log','XScale','log')

%% lat30 analysis

datstruct_lat30 = load("lat30.mat");
t = datstruct_lat30.time;
sum_insol = datstruct_lat30.summerinsol;

N = numel(t);
f = 0:1/N:1/2;

fft_si = fft(sum_insol);

Psi2 = abs(fft_si).^2/N;
Psi = Psi2(1:floor(N/2) +1);
Psi = [Psi(1);2*Psi(2:end)];

figure(14)
plot(f,Psi)
title('Summer Insolation PSD')
xlabel('Frequency, (cycles/kyr)')
set(gca,'YScale','log','XScale','log')