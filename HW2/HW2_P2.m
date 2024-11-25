close all
clear all
clc

%import data structure and extract dataspace
datstruct = load('boulder_temp.mat');
bulkdat = datstruct.temp;
% bulkdat column space names: Year, Month, Day, Temp (F).

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

figure(1)
tiledlayout(1,2)
%% ------------- Part A ------------------%%

%July Normal Probability Plot
nexttile
T_Jul = Tdata.Jul(:,4);
probplot(T_Jul)
title('July Normal Probability Plot')
xlabel('Temperature (F)')

%March Normal Probability Plot
nexttile
T_Mar = Tdata.Mar(:,4);
probplot(T_Mar)
title('March Normal Probaility Plot')
xlabel('Temperature (F)')

%% ------------ Part C ------------------------%

%Probability for temps less than 50 F in Nov, emperical and theoretical
figure(2)
T_Nov = Tdata.Nov(:,4);
%Emperical Distribution
[Nov_empT,Nov_edge,Nov_bin] = histcounts(T_Nov,'Normalization','probability','BinEdges',min(T_Nov):1:max(T_Nov));
% element 1 in Nov_empT is 15 F.
histogram(T_Nov,'Normalization','probability','BinEdges',min(T_Nov):1:max(T_Nov))
%Emperical probability
p_50emp = sum(Nov_empT(1:(49-14)));

%Theoretical Normal Distribution
x = min(T_Nov):1:max(T_Nov);
mu = mean(T_Nov);
s = std(T_Nov);
Nov_tT = normpdf(x,mu,s);
hold on
plot(x,normpdf(x,mu,s))
hold off
%Theoretical Probability
p_50t = sum(Nov_tT(1:(49-14)));

%Display answers
prompt1 = 'The probability of a temperature less than 50 F:';
tprompt1 = 'Theoretical =';
empprompt1 = 'Emperical =';

disp(prompt1)
disp(strcat(tprompt1,string(p_50t*100),'%'))
disp(strcat(empprompt1,string(p_50emp*100),'%'))

%% -------- Part D ---------- %%
mu_jul = mean(T_Jul);
sdevJul = std(T_Jul);
N_J = numel(T_Jul);
t_Jul = tinv(0.95,(N_J-1));
T_l = mu_jul-t_Jul*sdevJul/sqrt(N_J);
T_u = mu_jul+t_Jul*sdevJul/sqrt(N_J);

a = (T_l < T_Jul)&(T_Jul <T_u);
partd1 = sum(a)/N_J*100;

%get subset of july temperatures
ll = Tdata.Jul(:,1) == 1991 | Tdata.Jul(:,1) == 1992;
T_Jul92 = Tdata.Jul(ll,4);
mu_jul92 = mean(T_Jul92);
sdevJul92 = std(T_Jul92);
N_J92 = numel(T_Jul92);
t_jul92 = tinv(0.95,N_J92 - 1);
T_l92 = mu_jul92-t_jul92*sdevJul92/sqrt(N_J92);
T_u92 = mu_jul92+t_jul92*sdevJul92/sqrt(N_J92);
a92 = (T_l92 < T_Jul92)&(T_Jul92 <T_u92);
partd2 = sum(a92)/N_J92*100;


%% -----------Part E---------------%%
disp('Part E ---------------------------------------------------')
jj = Tdata.Jul(:,1) == 1997 | Tdata.Jul(:,1) == 1998 | Tdata.Jul(:,1) == 1999;
T_Jul9799 = Tdata.Jul(jj,4);

[h9799,p9799] = ztest(T_Jul9799,mu_jul,sdevJul);

if h9799 == 1
    disp('Tjul9799 is not a random sample')
else 
    disp('Tjul9799 is a random sample')
    disp(string(p9799))
end

kk = Tdata.Jul(:,1) == 2006 | Tdata.Jul(:,1) == 2007 | Tdata.Jul(:,1) == 2008; % adjust this line for [06,08]
T_Jul0607 = Tdata.Jul(kk,4);
[h0607,p0607] = ztest(T_Jul0607,mu_jul,sdevJul);

if h0607 == 1
    disp('Tjul0607 is not a random sample')
    
else
    disp('Tjul0607 is a random sample')
    disp(string(p0607))
end

%% ---------------- Part F --------------- %%
disp('Part F---------------------------------------')
T_Aug = Tdata.Aug(:,4);
mu_Aug = mean(T_Aug);
ttestaj = ttest2(T_Aug,T_Jul,"Alpha",0.05);
if ttestaj == 0
    disp('Both August and July have equal distributions')
else 
    disp('August and July do not have the same distributions')
end
nj = numel(T_Jul);
na = numel(T_Aug);
nu_aj = na + nj -2;
sdevAug = std(T_Aug);
t_AJ = (mu_jul - mu_Aug) / sqrt(  (nj+na)/(nj*na) * ( (nj-1)*(sdevJul^2) + (na-1) * (sdevAug^2))/nu_aj  );


%% ---------------Part G --------------------- %%

disp('Part G -----------------------------------------')
nn = Tdata.Aug(:,1) == 1991 | Tdata.Aug(:,1) == 1992;
T_aug91 = T_Aug(nn);
nu_aj91 = numel(T_aug91) + numel(T_Jul92) -2;
mu_aug91 = mean(T_aug91);
ttestaj92 = ttest2(T_aug91,T_Jul92,"Alpha",0.09);

if ttestaj92 == 0
    disp('The shortened series of July and August represent the same distribution')
end


%% ---------UDF's------------- %%

% function for data grouping by month
function m = monsort(mon_num,bulkdat,column)
    b = bulkdat(:,column) == mon_num;
    m = bulkdat(b,:);
end