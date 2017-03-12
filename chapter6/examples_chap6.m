%% this provides some examples for chapter 6.
% please contact Geliang Zhang, through zhanggel@gmail.com if you got any
% question.
%% example 1: use particle filter to filter the order flow intensity.
clear all;
clc;
load 'resampled_order_flow_example_one_single_day.mat';  

% Here we loaded one example of resampled order flow on one day. The time
% units is 10 seonds. 
% The data loaded here is processed by the author from one day of S&P
% 500 E-mini futures tick data (13/Jan/2011). The orignal dataset is too large to be updated and also is restricted by the data provider for redistribution.   
% The knowledge of processing of such dataset and reconstructing the LOB is beyond
% the scope of this thesis (could be found from other literature and references from the data provider), \
% and is not related with the academic content of this thesis. 
% Therefore, we omit the pre-processing steps of the dataset for the clarity. 

hit_level =5; % use the first five levels in particle filter algorithm. 

order_flow(:,1) = sum(params(:,1:1+hit_level-1),2); % ASK flow increase in sum; 
order_flow(:,2) = abs(sum(params(:,11:(11+hit_level-1)),2) - params(:,41)); % ASK flow cancelled decrease in sum; 
order_flow(:,3) = sum(params(:,21:(21+hit_level-1)),2); % BID flow increase in sum; 
order_flow(:,4) = abs(sum(params(:,31:(31+hit_level-1)),2) - params(:,42)); % BID flow cancelled decrease in sum; 
order_flow(:,5) = params(:,41);  % Market buy order;
order_flow(:,6) = params(:,42);  % Market sell order; 


order_flow(order_flow<=0) = 5;    

[trading_signal(1,:),err] = func_pf_tvPossion_opt_special(abs(order_flow(:,1)),100); % this function calls the partile filter to filter the order flow as well as giving the err bounds. 

% it is worthwhile to note that the err bounds given here is relative
% tight. This suggests a potential improvement over existing parameters of
% the PF approach. 

figure();
plot(5700:6000, order_flow(5700:6000,1));
xlabel('\bf Index of time / interval: 10 seconds'); ylabel('\bf Unfiltered intensity for order flow 1');
legend('Unfiltered intensity of order flow','location','best');

figure();
plot(5700:6000,trading_signal(1,5700:6000),'linewidth',1.5); grid();
hold on;
plot(5700:6000,trading_signal(1,5700:6000)+err(1,5700:6000),'r*'); hold on;
plot(5700:6000,trading_signal(1,5700:6000)-err(1,5700:6000),'r*');
xlabel('\bf Index of time / interval: 10 seconds'); ylabel('\bf Filtered intensity for order flow 1');
axis([5700,6000, min((trading_signal(1,5700:6000))), max((trading_signal(1,5700:6000)))]);
legend('Filtered intensity of order flow','Standard Error','location','best');



%% example 2: use the PF and the limit order book information to provide predictions. 

% Here we loaded one example of resampled order flow and the resample limit order book on one day. The time
% units is 10 seonds. 

% The data loaded here is processed by the author from one day of S&P
% 500 E-mini futures tick data (13/Jan/2011). The orignal dataset is too large to be updated and also is restricted by the data provider for redistribution.   
% The knowledge of processing of such dataset and reconstructing the LOB is beyond
% the scope of this thesis (could be found from other literature and references from the data provider), \
% and is not related with the academic content of this thesis. 
% Therefore, we omit the pre-processing steps of the dataset for the clarity. 
clear all; clc;

load 'resampled_order_flow_example_one_single_day.mat';  
load 'resampled_limit_order_book_example_one_single_day.mat';  
hit_level =5; % use the first five levels in particle filter algorithm. 




order_flow(:,1) = sum(params(:,1:1+hit_level-1),2); % ASK flow increase in sum of top hit_levels.; 
order_flow(:,2) = abs(sum(params(:,11:(11+hit_level-1)),2) - params(:,41)); % ASK flow cancelled decrease in sum of top hit_levels. 
order_flow(:,3) = sum(params(:,21:(21+hit_level-1)),2); % BID flow increase in sum of top hit_levels. 
order_flow(:,4) = abs(sum(params(:,31:(31+hit_level-1)),2) - params(:,42)); % BID flow cancelled decrease in sum of top hit_levels.
order_flow(:,5) = params(:,41);  % Market buy order;
order_flow(:,6) = params(:,42);  % Market sell order; 

% these data cleaning steps below are to make sure the stability of the algorithms used for PF and the conditional probability calculation. 
lob(isnan(lob))=0;  
order_flow(order_flow<=0) = 5;    

for i = 1:6
    trading_signal(i,:) = func_pf_tvPossion_opt(abs(order_flow(:,i)),100);
end    

para_pf = abs([trading_signal(1,:)',trading_signal(2,:)'+trading_signal(5,:)',trading_signal(3,:)',trading_signal(4,:)'+trading_signal(6,:)']);
[m,n] = size(para_pf);
para_pf(:,2) = max(para_pf(:,1:2),[],2)+10*ones(m,1);
para_pf(:,4) = max(para_pf(:,3:4),[],2)+10*ones(m,1);

ask_v = sum(lob(:,3:2:(1+hit_level*2)),2);
bid_v = sum(lob(:,4:2:(2+hit_level*2)),2);
p_inc = nan(m,1);
for i = 1 : m 
    if mod(i,10)==0
        i./m*100, % shows the percentage of progress. 
    end
    p_inc(i,1) = alg_pdf_fftv7(bid_v(i),ask_v(i),para_pf(i,:));  % this is the obtained conditional price-increase probability over time, which can be used for various testing purposes.
end




