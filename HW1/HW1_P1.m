clear all
close all
clc

%import data structure and extract dataspace
datstruct = load('boulder_precip.mat');
bulkdat = datstruct.precip;
% bulkdat column space names: Year, Month, Day, Precip Amt.

%initialize month str vec
m_str = string(zeros(12,1));

%separate data into months (column 2) and assign to monthly structure, name : data.month (ex: data.Jan)
for ii = 1:12
    date_form = datetime(100,ii,10);
    m_str(ii) = string(month(date_form,'shortname'));
    Pdata.(m_str(ii)) = monsort(ii,bulkdat);
end


% calculate average monthly precip and std. dev
averages = zeros(12,1);
std_devs = zeros(12,1);
for ii = 1:12
    precip_vals = Pdata.(m_str(ii))(:,4);
    averages(ii) = avg(precip_vals);
    std_devs(ii) = stdev(precip_vals,averages(ii));
end
x = 1:12;
Month = m_str;
Average_Precip = averages;
Ptable = table(Month,Average_Precip,std_devs);
uitable(uifigure,"Data",Ptable{:,:},'ColumnName', ...
    Ptable.Properties.VariableNames);
%create plot layout
figure(1)
%set(gcf, 'Position', get(0, 'Screensize'));
tiledlayout(2,2)

%plot averages
nexttile
plot(x,averages)
plt_attrib(' Monthly Precip', 'Average', 'inches', m_str, 'Month')

%plot std. devs
nexttile
plot(x,std_devs)
plt_attrib(' Monthly Precip', 'Standard Deviation', 'inches', m_str, 'Month')

%construct matrix of precip in june, column space = years
yr_span = max(Pdata.Jun(:,1)) - min(Pdata.Jun(:,1));
jun_mat = zeros(30,yr_span);
yr_0 = Pdata.Jun(1,1);
for jj = 1:yr_span
    yr = yr_0 - 1 + jj;
    b = Pdata.Jun(:,1) == yr;
    nz = find(b);
    jun_mat(:,jj) = Pdata.Jun(nz,4);
end

% make rain occurence matrix
Jun_rainday = jun_mat > 0;

% create number of rainy days each yr in june vector
drj = zeros(yr_span,1);
for jj = 1:yr_span
    drj(jj) = sum(Jun_rainday(:,jj));
end

%create histogram
nexttile
histogram(drj, 'BinEdges',0:1:30)
plt_attrib(' Rainy Days in June', 'Number',0,0, 'Number of Days in June' )
%set(gca,xtick = 0:2:30, xticklabels = 0:2:30)

nexttile
%create cumulative frequency plot of # of rainy day occurences
histogram(drj)
counts = histcounts(drj);

cf_drj = cf(counts);
plot(cf_drj)
title('CF of Rainy Days in June')
xlabel('Days')
ylabel('Probability')

partE = strcat('Probability of Exactly 10 Days of Rain in June = ', ...
    string(100*(cf_drj(11)-cf_drj(10))), '%');
disp(partE)

partF = strcat('Probability of more than 10 Days of Rain in June =', ...
    string(100*(cf_drj(numel(cf_drj))-cf_drj(11))), '%');
disp(partF)


%%_________UDFs__________%%

% function for data grouping by month
function m = monsort(mon_num,bulkdat)
    b = bulkdat(:,2) == mon_num;
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

% plot attributes
function plt_attrib(str_obj,str_title,y_units_str,mon_str,x_label)
    if class(mon_str) == 'string'
        set(gca, 'XTick',1:1:12, 'XTickLabel',mon_str)
    end
    hold on
    title(strcat(str_title, ' of ', str_obj))
    xlabel(x_label)
    if y_units_str ~= 0 
        ylabel(strcat(str_title,' (',y_units_str,')'))
    else 
        ylabel(str_title)
    end
    hold off
end