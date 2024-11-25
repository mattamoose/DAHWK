clear
close all
clc

%import data structure and extract dataspace
datstruct_b= load('boulder.mat');
bldr_droptime = datstruct_b.droptime;
bldr_time = datstruct_b.t;
bldr_y = datstruct_b.y;

datstruct_p = load('palmer.mat');
plmr_droptime = datstruct_p.droptime;
plmr_time = datstruct_p.t;
plmr_y = datstruct_p.y;
% bulkdat row space names: drop time (datestr), y(m) positive is down,time (s).

%% Part A ------------------------------------------------------

% Create Vandermonde Matrix for Boulder data:
% Calculate x via \
A_b = zeros(numel(bldr_y(:,1)),3);
g_b = zeros(numel(bldr_y(1,:)),1);
y0_b = zeros(size(g_b));
v0_b = zeros(size(g_b));
for jj = 1:numel(g_b)
    for ii = 1:numel(bldr_y(:,1))
        A_b(ii,:) = [1 bldr_time(ii,jj) 0.5*bldr_time(ii,jj)^2]; 
    end
x = A_b\bldr_y(:,jj);
y0_b(jj) = x(1);
v0_b(jj) = x(2);
g_b(jj) = x(3);
end

% Create Vandermonde Matrix for Palmer data:
% Calculate x via \
A_p = zeros(numel(plmr_y(:,1)),3);
g_p = zeros(numel(plmr_y(1,:)),1);
y0_p = zeros(size(g_p));
v0_p = zeros(size(g_p));
for jj = 1:numel(g_p)
    for ii = 1:numel(plmr_y(:,1))
        A_p(ii,:) = [1 plmr_time(ii,jj) 0.5*plmr_time(ii,jj)^2]; 
    end
x = A_p\plmr_y(:,jj);
y0_b(jj) = x(1);
v0_b(jj) = x(2);
g_p(jj) = x(3);
end

%% Part B ---------------------------------------------------
%average g 
g_bavg = mean(g_b);
g_pavg = mean(g_p);

disp(strcat('Average g in Boulder= ',string(g_bavg), ' m/s^2'))
disp(strcat('Average g in Palmer= ',string(g_pavg), ' m/s^2'))

%% Part C ----------------------------------------------------
% average y0
y_bavg = mean(y0_b);
y_pavg = mean(y0_p);

% average v0
v_bavg = mean(v0_b);
v_pavg = mean(v0_p);

f_b = @(t) y_bavg + v_bavg.*t + 0.5*g_bavg.*t.^2;
f_p = @(t) y_pavg + v_pavg.*t + 0.5.*g_pavg.*t.^2;

t = 0:0.01:bldr_time(numel(bldr_time(:,1)),1);
fb = f_b(t);
fp = f_p(t);

figure(1)
plot(t,fb)
hold on
plot(t,fp)
title('Free Fall')
legend('Boulder',"Palmer")
xlabel('time (s)')
ylabel('displacement (m)')

%% Problem 2 --------------------------------------------------------
g = [g_b(1:128,:) g_p];
p = anova1(g);
disp(strcat('Results from anova1 show that p=',string(p)))

%% Problem 3 ---------------------------------------------------------
c_b = zeros(1,3,numel(bldr_y(1,:)));
c_p = zeros(1,3,numel(plmr_y(1,:)));
for ii = 1:numel(bldr_y(1,:))
    c_b(:,:,ii) = polyfit(bldr_time,bldr_y,2);
end
c_b(1,1,:)=2*c_b(1,1,:);

for jj = 1:numel(plmr_y(1,:))
    c_p(:,:,jj) = polyfit(plmr_time,plmr_y,2);
end
c_p(1,1,:) = 2*c_p(1,1,:);

gply_b = mean(c_b(1,1,:));
gply_p = mean(c_p(1,1,:));

Gdiff_b = abs(gply_b-g_bavg);
Gdiff_p = abs(gply_p-g_pavg);
disp(strcat('The value of g in Boulder using Polyfit =',string(gply_b)))
disp(strcat(['The magnitude of the difference between the polyfit and \' ...
    'values of g in boulder=',string(Gdiff_b)]))
disp(strcat('The value of g in Palmer using Polyfit =', string(gply_p)))
disp(strcat(['The magnitude of the difference between the polyfit and \' ...
    'values of g in Palmer =',string(Gdiff_p)]))

%% Problem 4 ----------------------------------------------------------
g = @(x) 9.78*(1 + 0.0053*sin(x)^2);

latB = deg2rad(40.1308);
latP = deg2rad(61.5972);

gt_b = g(latB);
gt_p = g(latP);

gtdiffb = gt_b - g_bavg;
gtdiffp = gt_p - g_pavg;

disp(strcat('Theoretical g in Boulder=', string(gt_b)));
disp(strcat('The difference from theoretical in Boulder =',string(gtdiffb)));
disp(strcat('Theoretical g in Palmer =',string(gt_p)));
disp(strcat('The difference from theoretical in Palmer =',string(gtdiffp)));

%% Problem 5 ------------------------------------------------------------


y_b = zeros(size(bldr_y));
y_p = zeros(size(plmr_y));
b_yavg = zeros(numel(bldr_y(:,1)),1);

for kk = 1:numel(bldr_y(1,:))
    b_yavg(kk) = mean(bldr_y(kk,:));
end

for ii = 1:numel(bldr_y(1,:))
    y_b(:,ii) = y0_b(ii) + v0_b(ii).*bldr_time(:,ii) + 0.5*g_b(ii).*bldr_time(:,ii).^2;
    
end
SSE_b = zeros(numel(bldr_y(1,:)),1);
SSR_b = zeros(numel(bldr_y(1,:)),1);
YE_b = (bldr_y - y_b).^2;
YR_b = (y_b - b_yavg).^2;
for ii = 1:numel(bldr_y(1,:))
    SSE_b(ii) = sum(YE_b(:,ii));
    SSR_b(ii) = sum(YR_b(:,ii));
end
R_b = SSR_b./(SSR_b + SSE_b);
nu_b = numel(bldr_y(:,1)) - 3;
RMSE_b = sqrt(SSE_b./nu_b);

avgSSEb = mean(SSE_b);
avgSSRb = mean(SSR_b);
avgRsqb = mean(R_b);
avgRMSEb = mean(RMSE_b);


% For Palmer

y_p = zeros(size(plmr_y));
p_yavg = zeros(numel(plmr_y(:,1)),1);

for kk = 1:numel(plmr_y(1,:))
    p_yavg(kk) = mean(plmr_y(kk,:));
end

for jj = 1:numel(plmr_y(1,:))
    y_p(:,jj) = y0_p(jj) + v0_p(jj).*plmr_time(:,jj) + 0.5*g_p(jj).*plmr_time(:,jj).^2;
end

SSE_p = zeros(numel(plmr_y(1,:)),1);
SSR_p = zeros(numel(plmr_y(1,:)),1);
YE_p = (plmr_y - y_p).^2;
YR_p = (y_p - p_yavg).^2;
for ii = 1:numel(plmr_y(1,:))
    SSE_p(ii) = sum(YE_p(:,ii));
    SSR_p(ii) = sum(YR_p(:,ii));
end
R_p = SSR_p./(SSR_p + SSE_p);
nu_p = numel(plmr_y(:,1)) - 3;
RMSE_p = sqrt(SSE_p./nu_p);

avgSSEp = mean(SSE_p);
avgSSRp = mean(SSR_p);
avgRsqp = mean(R_p);
avgRMSEp = mean(RMSE_p);

disp('Average fit statistics for Boulder:')
disp(strcat('SSE_B =',string(avgSSEb)))
disp(strcat('R_sqr_B =',string(avgRsqb)))
disp(strcat('RMSE_B =', string(avgRMSEb)))

disp('Average fit statistics for Palmer:')
disp(strcat('SSE_P =',string(avgSSEp)))
disp(strcat('R_sqr_P =',string(avgRsqp)))
disp(strcat('RMSE_P =', string(avgRMSEp)))

%% Problem 6 ------------------------------------------------------------
res_b = bldr_y - y_b;
nd_b = ones(numel(bldr_y(1,:)),1);

for i = 1:numel(bldr_y(1,:))
    m_rb = mean(res_b(:,i));
    s_rb = std(res_b(:,i));
    %nd_b(i) = ztest(res_b(:,i),m_rb,s_rb);
    %nd_b(i) = ttest(res_b(:,i),m_rb);
    nd_b(i) = chi2gof(res_b(:,i));
end

not_nd_b = find(nd_b);

res_p = plmr_y - y_p;
nd_p = ones(numel(plmr_y(1,:)),1);

for i = 1:numel(plmr_y(1,:))
    x = res_p(:,i);
    m_rp = mean(x);
    s_rp = std(x);
    %nd_p(i) = ztest(res_p(:,i),m_rp,s_rp);
    %nd_p(i) = ttest(res_p(:,i),m_rp);
    nd_p(i) = chi2gof(x);
end

not_nd_p = find(nd_p);
disp('Problem 6')
disp('For Boulder')
disp(not_nd_b')
disp('for Palmer')
disp(not_nd_p')
%% Problem 7 -------------------------------------------------------------

res_bcom = reshape(res_b,[],1);
res_pcom = reshape(res_p,[],1);

[h_b, p_b] = chi2gof(res_bcom);
[h_p, p_p] = chi2gof(res_pcom);

disp('Problem 7')
disp('For Boulder:')
disp(h_b)
disp(p_b)
disp('For Palmer:')
disp(h_p)
disp(p_p)

figure(4)
tiledlayout(1,2)
nexttile
histogram(res_bcom)
title('Boulder Regression Residual Distribution')
xlabel('Residual')
ylabel('Frequency')
nexttile
probplot(res_bcom)
title("Probability Plot of Boulder Regression Residuals")

figure(5)
tiledlayout(1,2)
nexttile
histogram(res_pcom)
title('Palmer Regression Residual Distribution')
xlabel('Residual')
ylabel('Frequency')
nexttile
probplot(res_pcom)
title('Probability Plot of Palmer Regression Residuals')

%% Problem 8 -------------------------------------------------------------

% Create Vandermonde Matrix for Boulder data:
% Calculate x via QR
A_b = zeros(numel(bldr_y(:,1)),3);
g_bqr = zeros(numel(bldr_y(1,:)),1);
y0_bqr = zeros(size(g_b));
v0_bqr = zeros(size(g_b));
for jj = 1:numel(g_b)
    for ii = 1:numel(bldr_y(:,1))
        A_b(ii,:) = [1 bldr_time(ii,jj) 0.5*bldr_time(ii,jj)^2]; 
    end
[Q,R] = qr(A_b,0);
x = R\(Q'*bldr_y(:,jj));
y0_bqr(jj) = x(1);
v0_bqr(jj) = x(2);
g_bqr(jj) = x(3);
end

% Create Vandermonde Matrix for Palmer data:
% Calculate x via QR
A_p = zeros(numel(plmr_y(:,1)),3);
g_pqr = zeros(numel(plmr_y(1,:)),1);
y0_pqr = zeros(size(g_p));
v0_pqr = zeros(size(g_p));
for jj = 1:numel(g_p)
    for ii = 1:numel(plmr_y(:,1))
        A_p(ii,:) = [1 plmr_time(ii,jj) 0.5*plmr_time(ii,jj)^2]; 
    end
[Q,R] = qr(A_p,0);
x = R\(Q'*plmr_y(:,jj));
y0_bqr(jj) = x(1);
v0_bqr(jj) = x(2);
g_pqr(jj) = x(3);
end

g_bqravg = mean(g_bqr);
g_pqravg = mean(g_pqr);

disp(strcat('Using QR, average g for Boulder = ',string(g_bqravg)))
disp(strcat('Using QR, average g for Palmer =',string(g_pqravg)))

%% Problem 9 --------------------------------------------------------------

% Create Vandermonde Matrix for Boulder data:
% Calculate x via QR
A_b = zeros(numel(bldr_y(:,1)),3);
g_bsvd = zeros(numel(bldr_y(1,:)),1);
y0_bsvd = zeros(size(g_b));
v0_bsvd = zeros(size(g_b));
for jj = 1:numel(g_b)
    for ii = 1:numel(bldr_y(:,1))
        A_b(ii,:) = [1 bldr_time(ii,jj) 0.5*bldr_time(ii,jj)^2]; 
    end
[U,Sb,V] = svd(A_b,0);
sr = zeros(size(Sb'));
St = Sb';
for kk = 1:3
    sr(kk,kk) = 1/St(kk,kk);
end
sAinv = V*sr*U';
x = sAinv*bldr_y(:,jj);
y0_bsvd(jj) = x(1);
v0_bsvd(jj) = x(2);
g_bsvd(jj) = x(3);
end

% Create Vandermonde Matrix for Palmer data:
% Calculate x via QR
A_p = zeros(numel(plmr_y(:,1)),3);
g_psvd = zeros(numel(plmr_y(1,:)),1);
y0_psvd = zeros(size(g_p));
v0_psvd = zeros(size(g_p));
for jj = 1:numel(g_p)
    for ii = 1:numel(plmr_y(:,1))
        A_p(ii,:) = [1 plmr_time(ii,jj) 0.5*plmr_time(ii,jj)^2]; 
    end
[U,S,V] = svd(A_p,0);
sr = zeros(size(S'));
St = S';
for kk = 1:3
    sr(kk,kk) = 1/St(kk,kk);
end
sAinv = V*sr*U';
x = sAinv*plmr_y(:,jj);
y0_bsvd(jj) = x(1);
v0_bsvd(jj) = x(2);
g_psvd(jj) = x(3);
end

g_bsvdavg = mean(g_bsvd);
g_psvdavg = mean(g_psvd);

disp(strcat('Using SVD, average g for Boulder = ',string(g_bsvdavg)))
disp(strcat('Using SVD, average g for Palmer =',string(g_psvdavg)))