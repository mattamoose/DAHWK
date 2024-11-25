close all; clear all; clc

bulkdat = load('globalsl.mat');
data.y = bulkdat.sealevel;
data.t = bulkdat.time;

figure(1)

scatter(data.t,data.y)
title('Sea Level vs Time')
xlabel('time (years)')
ylabel('Sea Level')


mdl = fitlm(data.t,data.y)
b = mdl.Coefficients.Estimate;
SE = mdl.Coefficients.SE


y_lm = @(x) b(2)*x + b(1);


hold on
plot(data.t, y_lm(data.t))

CI = coefCI(mdl);

disp('CI')
disp((CI))

resid = data.y - y_lm(data.t);

figure(2)
plot(data.t,resid)
title('Regression Residuals vs time')
figure(4)
histogram(resid)

[h1,p1] = chi2gof(resid)

%% Part e
n = numel(data.t);
H1 = ones(n,4);


H1(:,2:4) = [data.t cos(2*pi*data.t) sin(2*pi*data.t)];


be = H1 \ data.y

figure(3)
scatter(data.t,data.y)
title('Sea Level vs Time')
xlabel('time (years)')
ylabel('Sea Level')
hold on
plot(data.t,H1*be,LineWidth=2)


H2 = ones(n,6);

H2(:,2:6) = [data.t cos(2*pi*data.t) sin(2*pi*data.t) cos(4*pi*data.t) sin(4*pi*data.t)];

bf = H2 \data.y
plot(data.t,H2*bf,LineWidth=2)
legend('','Single Mode','Dual Mode')