

%import data structure and extract dataspace
datstruct = load('boulder_temp.mat');
bulkdat = datstruct.temp;
% bulkdat column space names: Year, Month, Day, Temp.

%initialize month and yr str vec
m_str = string(zeros(12,1));
yr_0 = bulkdat(1,1);
yr_f = bulkdat(numel(bulkdat(:,1)),1);
yr_str = string(zeros((yr_f-yr_0),1));
%separate data into years (column 1) and months (column 2) and assign to monthly structure, name : data.month (ex: data.Jan)
for ii = 1:12
    date_form = datetime(100,ii,10);
    m_str(ii) = string(month(date_form,'shortname'));
    Tdata.(m_str(ii)) = monsort(ii,bulkdat,2);
end

for jj = 0:(yr_f - yr_0)
    yr = yr_0 + jj;
    yr_str(jj + 1) = string(yr);
    fm = strcat('yr',string(yr));
    Tdata.(fm) = monsort(yr,bulkdat,1);
end



% get yearly averages
yravg = zeros(numel(yr_str),1);

for kk = 1:numel(yravg)
    fm = strcat('yr',yr_str(kk));
    yravg(kk) = avg(Tdata.(fm)(:,4));
end

% project each yearly average onto whole year in bulkdat
avgcol = zeros(numel(bulkdat(:,1)),1);
for ll = 1:numel(yr_str)
    a = string(bulkdat(:,1)) == yr_str(ll);
    avgcol(a) = yravg(ll);
end


% Temperature Vector
T = bulkdat(:,4);
days = 1:numel(bulkdat(:,1));

%create plot layout
figure(2)
%set(gcf, 'Position', get(0, 'Screensize'));
tiledlayout(2,2)

% plot Daily Temperatures
nexttile
plot(days,T)
hold on
plot(days,avgcol)
hold off
legend('Daily Temperature','Yearly Average')
xticks(0:500:numel(T))
xlabel('Time (Days)')
ylabel('Temperature (F)')
title('Temperatures since 1991 Jan 1')


% Plot July Data in Histogram
Tjuly = Tdata.Jul(:,4);
nexttile
histogram(Tjuly)
title('July Temperature Distribution')
ylabel('Number of Days')
xlabel('Days of July')
%CFD of july
counts = histcounts(Tjuly);
Tcfd = cf(counts);

%plot CFD of July
nexttile
Tjmin = min(Tjuly);
Tjmax = max(Tjuly);
x = Tjmin : Tjmax;

plot(x,Tcfd)
p_80 = Tcfd(x == 80);
title('CFD of July Temperatures')
xlabel('Temperature (F)')
ylabel('Probability')
answer = strcat(['Probability of temperature in July greater than ' ...
    '80 F = '],string( 100*(1-p_80)),'%');
disp(answer);

%march data
march_temp = Tdata.Mar(:,4);
nexttile
histogram(march_temp,'BinEdges',[min(march_temp):1:max(march_temp)])
%xticks(min(march_temp):5:max(march_temp))
title('Temperature Distribution in March')
xlabel('Temperature (F)')
ylabel('Number of Days')

skew = skewness(march_temp)
kurt = kurtosis(march_temp)

maxT = max(T)
minT = min(T)

%%_________UDFs__________%%

% function for data grouping by month
function m = monsort(mon_num,bulkdat,column)
    b = bulkdat(:,column) == mon_num;
    m = bulkdat(b,:);
end

% function to calculate mean
function mean = avg(vec)
    tot = sum(vec);
    mean = tot / numel(vec);
end

% function to calculate std.dev
function s = stdev(mon_vec,mon_avg)
    s = sqrt(1 / (numel(mon_vec) - 1) * sum((mon_vec - mon_avg).^2));
end

%create cumulative frequency data
function c = cf(vec)
   tot = sum(vec);
   c = zeros(numel(vec),1);
   for i = 1:numel(vec)
        relsum = sum(vec(1:i));
        c(i) = relsum/tot;
   end
end