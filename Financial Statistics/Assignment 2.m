%% Plot hourly log-close price and log-return
lreturns = diff(EOS.logclose);
logcl = EOS.logclose;
figure
plot(logcl)
title('Log-close price vs Time (EOS)')
xlabel('Time (hours)')
ylabel('log-close price')
figure
plot(lreturns)
title('Log-return vs Time (EOS)')
xlabel('Time (hours)')
ylabel('log-return')

[a1,lags1] = autocorr(lreturns,250);
plot(lags1,a1,'-r')
hold on
[a2,lags2] = autocorr(abs(lreturns),250);
plot(lags2,a2,'-m')
axis([0 250 -.1 1])
xlabel('lags','fontsize',14)
ylabel('autocorrelation','fontsize',14)
title('Autocorrelation (EOS)','fontsize',250)
set(gca,'fontsize',14)
legend({'log-returns','|log-returns|'})

%% Aggregation of returns 

tau = 24; % Number of days for aggregation of returns 

flag = 1; % Flag variable (set to 1 to aggregate returns over longer time scales)

if flag == 1
    
    aux = [];
   
    for t = 0:tau:length(lreturns)-tau
       
        aux = [aux; sum(lreturns(t+1:t+tau))];
        
    end
    
    lreturns = aux;
    
end

plot(lreturns)
%% First four moments 

% Compute and print the values of the first four moments

N = length(lreturns);
m = sum(lreturns)/N; % Compute mean and store value in variable
fprintf('\n')
fprintf('Mean = %4.3f\n',m)
s = sqrt(sum((lreturns-m).^2)/N); % Compute std. deviation and store value in variable
fprintf('Std. deviation = %4.3f\n',s)
fprintf('Skewness = %4.3f\n',sum((lreturns-m).^3)/(N*s^3))
fprintf('Excess kurtosis = %4.3f\n',sum((lreturns-m).^4)/(N*s^4)-3)
fprintf('\n')

%% Plot relative frequency distribution of returns
% compute and plot return distribution
[freq,bin]=hist(lreturns,100);
figure
bar(bin,freq/sum(freq))
hold on
% compare with normal distribution
x = min(lreturns):(max(lreturns)-min(lreturns))/1000:max(lreturns);
plot(x,normpdf(x,m,s)*(bin(2)-bin(1)),'-m','linewidth',2)
axis([-0.2 0.2 0 max(freq/sum(freq))*1.1])
legend({'log-return','normal'})
title({'Log-return relative frequency distribution (EOS)'})
xlabel('log-return')
ylabel('relative frequency')

%% qq plot

qqplot(lreturns);
title('QQ Plot of EOS vs Standard Normal ')

%% Plot ccdf of returns 
loglog(sort(lreturns(lreturns>0)),1-[1:(length(lreturns(lreturns>0)))]/length(lreturns(lreturns>0)),'+b');
hold on
loglog(sort(-lreturns(lreturns<0)),1-[1:(length(lreturns(lreturns<0)))]/length(lreturns(lreturns<0)),'xr');
% compare with normal distribution
x = max(lreturns)/1000:max(lreturns)/1000:max(lreturns);
loglog(x,1-(normcdf(x,m,s)-0.5)*2,'-m','linewidth',2)
axis([1e-5 0.5 1e-4 1])
legend('pos ret','neg ret','normal', 'Location', 'northwest')
title('Complemetary cumulative distribution (EOS)')
xlabel('log-return')
ylabel('complemetary cumulative distribution')

%% Fitting Right and Left Tail via Maximum Likelihood and Bootstrap

p = 0.1; % Defining tails as top p% of returns (both positive and negative)

%%% Right tail
%% GEV
lreturns = sort(lreturns); % Sorting returns
r_right = lreturns(round((1-p)*length(lreturns)):end); % Selecting top p% of returns
[F,r_righti]=ecdf(r_right);
parmhat = gevfit(r_right);

plot(r_righti,gevcdf(r_righti,0.6128,0.0072,0.0281),'-')
hold on
stairs(r_righti,F,'r')
legend('Fitted GEV CDF', 'Empirical CDF','Location','best')


test_cdf = makedist('GeneralizedExtremeValue','k',0.6128,'sigma',0.0072,'mu',0.0281);
[h,p] = kstest(r_right,'CDF',test_cdf);

%% Left tail

r_left = lreturns(1:round(p*length(lreturns))); % Selecting bottom p% of returns
r_left = abs(r_left); % Converting negative returns to positive numbers
[F,r_lefti]=ecdf(r_left);
parmhat = gevfit(r_left);

plot(r_lefti,gevcdf(r_lefti,0.6371,0.0063,0.0269),'-')
hold on
stairs(r_lefti,F,'r')
legend('Fitted GEV CDF', 'Empirical CDF','Location','best')
title('EOS & GEV comparison for tail index')



test_cdf = makedist('GeneralizedExtremeValue','k',0.6371,'sigma',0.0063,'mu',0.0269);
[h,p] = kstest(r_left,'CDF',test_cdf);

%% Kernel density
[f,xi,bw] = ksdensity(lreturns,'Bandwidth',0.03); 
bw
figure
plot(xi,f,'--r','LineWidth',1.0);
title({'Kernel density estimation (EOS)';'(bandwidth = 0.03)'})
xlabel('log-return')


[f,xi,bw] = ksdensity(lreturns,'Bandwidth',0.03,'Function','cdf'); %'Bandwidth',0.3,); 
bw
figure
plot(xi,f,'--b','LineWidth',1.0);
title({'Cumulative kernel density (EOS)';'(bandwidth = 0.03)'})
xlabel('log-return')

%% VaR,cVaR and bootstrap
VaR = quantile(lreturns,0.05);
cVaR = mean(lreturns(lreturns<=VaR));

%r = sort(lreturns); % Sorting returns
Nbts = 1000; % Number of bootstrap samples
alpha = 0.95;
bts = 0.8; % Fraction of data to be retained in each bootstrap sample

 bts_VaR = []; 
 bts_cVaR = [];
 
for i = 1:Nbts
    r_bts = lreturns(randperm(length(lreturns))); % Random permutation of returns
    r_bts = r_bts(1:round(bts*length(r_bts))); % Bootstrapping bts% of returns 
    VaR = quantile(r_bts,0.05);
    cVaR = mean(r_bts(r_bts<=VaR));
    
    bts_VaR = [bts_VaR; VaR];
    bts_cVaR = [bts_cVaR; cVaR];
end

bts_VaR = sort(bts_VaR); % Sorting bootstrap estimates for right tail exponent
bts_cVaR = sort(bts_cVaR); % Sorting bootstrap estimates for right tail exponent
fprintf('VaR interval at %3.2f CL: [%4.3f; %4.3f] \n',alpha, bts_VaR(round(0.5*(1-alpha)*Nbts)), bts_VaR(round(0.5*(1+alpha)*Nbts)))
fprintf('\n')
fprintf('cVaR interval at %3.2f CL: [%4.3f; %4.3f] \n',alpha, bts_cVaR(round(0.5*(1-alpha)*Nbts)), bts_cVaR(round(0.5*(1+alpha)*Nbts)))

%%Volume vs log return (EOS)

vol = EOS.total_volume_global
scatter(vol,logcl,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
title('Volume vs log-return (EOS)')
xlabel('Volume')
ylabel('log-return')
%% Plot both ccdf and power law tails

% r = polyfit(log(sort(-lreturns(lreturns<0))),log(1-[1:(length(lreturns(lreturns<0)))]/length(lreturns(lreturns<0)))); 
% %yApprox = (10^intercept)*x.^(slope);
% loglog(sort(-lreturns(lreturns<0)),1-[1:(length(lreturns(lreturns<0)))]/length(lreturns(lreturns<0)),'xr');

% [ycdf,xcdf] = cdfcalc(-lreturns(lreturns<0));
% xccdf = xcdf;
% yccdf = 1-ycdf(1:end-1);
% loglog(xccdf,yccdf,'xr');
% hold on
% 
% p = polyfit(log(xccdf),log(yccdf),1);
% m = p(1);
% b = exp(p(2));
% loglog(xccdf, b*xccdf.^m)
% 
% % y = r(2)*xccdf.^r(1);
% % plot(log(xccdf),log(y),'-b'); 
% %legend('data','linear fit')
% % % compare with normal distribution
% % x = max(lreturns)/1000:max(lreturns)/1000:max(lreturns);
% % %axis([1e-5 0.5 1e-4 1])
% % 
% % 
% % legend('neg ret','ae^{\alpha}', 'Location', 'best')
% % title('Complemetary cumulative distribution (Ripple)')
% % xlabel('log-return')
% % ylabel('complemetary cumulative distribution')
% r = sort(lreturns); % Sorting returns
% r_left = r(1:round(p*length(r))); % Selecting bottom p% of returns
% r_left = abs(r_left); % Converting negative returns to positive numbers
% 
% N = length(r_left); % Number of returns selected as left tail
% alpha_left = N/sum(log(r_left/min(r_left))); % Maximum-likelihood estimate for left tail exponent
% 
% fprintf('Left tail exponent: %4.3f\n',alpha_left)
% 
% x_left = linspace(min(r_left),max(r_left),100);
% y_left = alpha_left*(x_left/min(r_left)).^(-alpha_left-1)/min(r_left); % Power law distribution
% loglog(x_left,y_left,'r','LineWidth',2)
% set(gca,'FontSize',20)
% title('Left tail')
% 



