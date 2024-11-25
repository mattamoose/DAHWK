close all; clear all; clc

bulkdat = load('coloairquality.mat');

denver.time = bulkdat.airqual(:,2);
boulder.time = bulkdat.airqual(:,2);
longmont.time = bulkdat.airqual(:,2);

for i = 3:23
    fname = string(bulkdat.fnames1(i));
    if i >= 3 && i <= 9
        boulder.(fname) = bulkdat.airqual(:,i);
    elseif i >= 10 && i <= 16
        longmont.(fname) = bulkdat.airqual(:,i);
    elseif i >= 17 && i <= 23
        denver.(fname) = bulkdat.airqual(:,i);
    end
end

%% Part A

% get components of ws direction

% east-west => x => cos()
% north-west => y => sin()
wdrad = deg2rad(denver.WDdenver(:));
denver.ew = cos(wdrad);
denver.ns = sin(wdrad);
b0 = ones(numel(denver.ew(:,1)),1);
obs = [b0 ,denver.COdenver, denver.COxdenver, denver.WSdenver, denver.ew,denver.ns, denver.Tdenver];

obs1 = obs;
obs2 = obs;
[b1, bint1, r1, rint1,stats1] = regress(denver.O3denver, obs);
disp(b1)
disp(bint1)
ym = obs*b1 - denver.O3denver;
obs1(:,4) = [];
[b2, bint2,r2,rint2,stats2] = regress(denver.O3denver, obs1);

obs2(:,5:6) = [];
[b3, bint3,r3,rint3,stats3] = regress(denver.O3denver,obs2);
ym3 = obs2*b3 - denver.O3denver;
[h1,p1] = kstest2(r1.^2,r3.^2)

disp('The wind direction is statistically insignificant')



%% Part c

b0 = regress(denver.O3denver,obs);
yfit = obs*b0;
res = denver.O3denver - yfit;
ci = bootci(1000,{@(x)regress(yfit + x,obs),res},'type','normal')


%% Part D

%calculate correlation between adjacent measurements of ozone, then points
%24 hrs apart; for all 3 locations.

corrD = sqrt(stats1(1,1))
wbrad = deg2rad(boulder.WDboulder(:));
denver.ew = cos(wbrad);
denver.ns = sin(wbrad);
b0 = ones(numel(boulder.O3boulder),1);
obsB = [b0 ,boulder.COboulder, boulder.COxboulder, boulder.WSboulder, denver.ew,denver.ns, boulder.Tboulder];
wdrad = deg2rad(longmont.WDlongmont(:));
denver.ew = cos(wdrad);
denver.ns = sin(wdrad);
obsL = [b0, longmont.COlongmont, longmont.COxlongmont, longmont.WSlongmont, denver.ew, denver.ns, longmont.Tlongmont];

[bb, bintb, rb, rintb,statsb] = regress(boulder.O3boulder, obsB);
[bl, bintl, rl, rintl, statsl] = regress(longmont.O3longmont,obsL);

corrB = sqrt(statsb(1,1))
corrL = sqrt(statsl(1,1))

