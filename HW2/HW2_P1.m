 close all
clear all
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

%% April Portion of Problem 1

%construct matrix of precip in April, column space = years
yr_span = max(Pdata.Apr(:,1)) - min(Pdata.Apr(:,1));
apr_mat = zeros(30,yr_span);
yr_0 = Pdata.Apr(1,1);
for jj = 1:yr_span
    yr = yr_0 - 1 + jj;
    b = Pdata.Apr(:,1) == yr;
    nz = find(b);
    apr_mat(:,jj) = Pdata.Apr(nz,4);
end

% make rain occurence matrix
apr_rainday = apr_mat > 0;

% create number of rainy days each yr in Apr vector
dra = zeros(yr_span,1);
for jj = 1:yr_span
    dra(jj) = sum(apr_rainday(:,jj));
end


%create histogram for Apr
nexttile
histogram(dra,'Normalization','probability', 'BinEdges',0:1:30)
plt_attrib(' Rainy Days in April', 'Probability',0,0, 'Number of Days in April' )
[aprp_emp,edges,bin] = histcounts(dra,'Normalization','probability', 'BinEdges',0:1:30);


%calculate and plot the binomial distribution using UDF.
p = sum(dra)/numel(apr_mat);
n = 30;
bd = bi_dist(p,n);

hold on
plot(bd)

%Print answers for part D
disp('Based on the emperical probability distribution')
ans_prompt1 = 'The probability of 10 days of rain in April =';
ans_de = strcat(ans_prompt1, string(aprp_emp(11)*100),'%');
disp(ans_de)



disp('Based on the theoretical binomial probability distribution')
ans_dt = strcat(ans_prompt1, string(bd(11)*100), '%');
disp(ans_dt)


%% Part E

%construct matrix of precip in May, column space = years
yr_span = max(Pdata.May(:,1)) - min(Pdata.May(:,1));
may_mat = zeros(31,yr_span);
yr_0 = Pdata.May(1,1);
for jj = 1:yr_span
    yr = yr_0 - 1 + jj;
    b = Pdata.May(:,1) == yr;
    nz = find(b);
    may_mat(:,jj) = Pdata.May(nz,4);
end

% make rain occurence matrix
may_rainday = may_mat > 0;

% create number of rainy days each yr in May vector
drm = zeros(yr_span,1);
for jj = 1:yr_span
    drm(jj) = sum(may_rainday(:,jj));
end

% get emperical and theoretical probability for 10 days
mayp_emp = histcounts(drm,'Normalization','probability');
p_10emp = cumbidist(mayp_emp,11);

n_m = 31;
p_m = sum(drm)/numel(may_mat);
m_bd = bi_dist(p_m,n_m);
p_10t = cumbidist(m_bd,11);

% plot answers
disp('Based on the emperical cummulative probability')
ans_prompt1 = 'The probability of more than 10 days of rain in May =';
ans_de = strcat(ans_prompt1, string((1 - p_10emp)*100),'%');
disp(ans_de)



disp('Based on the theoretical binomial probability distribution')
ans_dt = strcat(ans_prompt1, string((1 - p_10t)*100), '%');
disp(ans_dt)


%%
%%%%%%%%%% User Defined Functions %%%%%%%%%%%%%%%

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

% function for data grouping by month
function m = monsort(mon_num,bulkdat)
    b = bulkdat(:,2) == mon_num;
    m = bulkdat(b,:);
end

% binomial distribution PDF
function bd = bi_dist(p,n)
    bd = zeros(n,1);
    for i = 1:30
        a = factorial(n)/(factorial(i)*(factorial(n-i)));
        bd(i) = a*p^i*(1-p)^(n-i);
    end
    
end

% binomial CDF
function cbd = cumbidist(d,i)
    v = resize(d,i);
    cbd = sum(v);
end