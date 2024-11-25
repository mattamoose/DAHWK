
data = load("homework8.mat");
t = data.time; lon = (data.lon); lat = double(data.lat); sl = double(data.sldata);
% keep = {'t','lon','lat','sl'};
% clear(setdiff(who,keep)); clc; close all
clc;close all
%spatial res 1x1 degree lat-lon units = cm

%sl(i,j,:) = [];

MO = NaN(76,360*180);
ll = 1;
for ii = 1:numel(lon)
    for jj = 1:numel(lat)
      k = sl(ii,jj,:);
      MO(:,ll) = k;
      ll = ll+1;
    end
end

land = 998*ones(76,1); A = MO; l = 0; island_col = NaN(ii*jj,1);
for kk = 1:ii*jj

    if A(:,kk) == land
        MO(:,kk-l) = [];
        island_col(kk) = kk;
        l = l+1;
    end
end

%% 
island_col(isnan(island_col)) = [];

% interp land
B = MO; % all steady state land removed
still_land = find(B == 998);

[land_row,land_col] = find(B == 998);

Nr = size(MO,1); Nc = size(MO,2);
x = 1:76; xq = x;
for jj = 1:Nc
    v = MO(:,jj);
    lnd = find( v == 998);
    x(lnd) = []; v(lnd) = [];
    
    MO(:,jj) = interp1(x,v,xq,'linear','extrap');
    x = 1:76;
end

check1 = find(isnan(MO)); check2 = find(MO == 998);
%disp(check1)
disp(check2)
%% Centering the data
for jj = 1:Nc
    v = MO(:,jj);
    MO(:,jj) = v - mean(v);

end
disp('data is centered')

%% SVD

[u,s,v] = svd(MO,'econ');

lambda = s.^2;
vtot = sum(lambda,'all');
lambda = 100*lambda/vtot;

disp(lambda(1,1))
disp(lambda(2,2))
disp(lambda(3,3))
disp(lambda(4,4))
disp(lambda(5,5))

V1 = v*s;

for i = 1:size(v,2)
    Vmean(i) = mean(V1(:,i));
end

for j = 1:size(u,2)
    PCTS(:,j) = u(:,j)*Vmean(j);
end

num_modes = 5;
LVm = zeros(size(v));

for k = 1:size(v,2)
    LVm(:,k) = V1(:,k)./Vmean(k);
end

%%

EOF = PCTS(:,5)*LVm(:,5)';
EOF5 = EOF;


for i = 1:num_modes
    titl = strcat('PCTS mode',string(i));
    figure(i)
    plot(t,PCTS(:,i))
    title(titl)
end

%% Add back land to LVm

Lvm1 = LVm(:,1); Lvm2 = LVm(:,2); Lvm3 = LVm(:,3);
Lvm4 = LVm(:,4); Lvm5 = LVm(:,5);


% Lvm1 = iso(land_col,land_row,Lvm1);
% Lvm2 = iso(land_col,land_row,Lvm2);
% Lvm3 = iso(land_col,land_row,Lvm3);
% Lvm4 = iso(land_col,land_row,Lvm4);
% Lvm5 = iso(land_col,land_row,Lvm5);


Lvm1 = island(island_col,Lvm1);
Lvm2 = island(island_col,Lvm2);
Lvm3 = island(island_col,Lvm3);
Lvm4 = island(island_col,Lvm4);
Lvm5 = island(island_col,Lvm5);

Lvm1 = res(Lvm1);
Lvm2 = res(Lvm2);
Lvm3 = res(Lvm3);
Lvm4 = res(Lvm4);
Lvm5 = res(Lvm5);


figure(6)
contour(Lvm1')
title('mode 1')

figure(7)
contour(Lvm2')
title('mode 2')

figure(8)
contour(Lvm3')
title('mode 3')

figure(9)
contour(Lvm4')
title('mode 4')

figure(10)
contour(Lvm5')
title('mode 5')

%% Periodogram of time series
figure(11)
periodogram(PCTS(:,1))
title('mode 1')

figure(12)
periodogram(PCTS(:,2))
title('mode 2')

figure(13)
periodogram(PCTS(:,3))
title('mode 3')
%% Add back in land

pca = EOF5;

island_col(isnan(island_col)) = [];

for ii = 1:numel(island_col)
    jj = island_col(ii);
    if jj == 1
        a = pca;
        pca = [land a];
    elseif jj == size(pca,2)
        pca = [pca land];
    else 
        a = pca(:,1:jj -1); b = pca(:,jj:end);
        pca = [a land b];
    end
end
%%
% kappa = pca;

%% 
clear alpha gamma
pca = kappa;

beta = NaN(360,180,76);
ll =1;mm = 0; jj = 1;
for ii = 1:size(pca,2)
    kk = ii - mm;
    beta(jj,kk,:) = pca(:,ii);
    if ll == 180
        mm = mm + 180;
        ll = 0;
        jj = jj + 1;
    end
    ll = ll +1;
    
end

alpha = beta(:,:,1);
gamma = sl(:,:,1);

diff = alpha - gamma;


%% recontructed visualization

figure(14)
tiledlayout (2,1)
nexttile
contourf(beta(:,:,76)'), colorbar
title('reconstructed data')

nexttile
contourf(sl(:,:,76)'), colorbar
title('original data')

%% Functions

function M = iso(land_col,land_row,Lv)
    for ii = 1:numel(land_col)
    disp(ii)
    j = land_col(ii); k = land_row(ii);
    i = sub2ind(size(Lv),k,j);
    Lv(i) = 998;
    
    end
    M = Lv;
end

function M = island(island_col,V)
    for ii = 1:numel(island_col)
        jj = island_col(ii);
        if jj == 1
            a = V;
            V = [998; a];
        elseif jj == size(V,1)
            V = [V; 998];
        else 
            a = V(1:jj -1); b = V(jj:end);
            V = [a; 998; b];
        end
    end
    M = V;
end

function M = res(V)
M = NaN(360,180);
ll =1;mm = 0; jj = 1;
for ii = 1:size(V,1)
    kk = ii - mm;
    M(jj,kk) = V(ii);
    if ll == 180
        mm = mm + 180;
        ll = 0;
        jj = jj + 1;
    end
    ll = ll +1;
    
end
end