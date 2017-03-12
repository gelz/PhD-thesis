function [returnData] = func_estimate_statistics()

maindir = 'C:\es\';

files_path =  fuf([maindir,'\*preRun.mat'],'detail');  % here to change the file paths for the E-mini futures dataset.

table = nan(length(files_path),41);
% for i = 1:length(files_path)
for i = 1:length(files_path)

    try
        returnData{i,1} = local_func_find_params(files_path{i});
    catch me
        continue;
    end
    
    names = fieldnames(returnData{i,1});
    
    for j = 1: length(names)
        table(i,j) = getfield(returnData{i,1},names{j});
    end
    
    
%     keyboard;
end

final_table(:,1) = nanmean(table,1);
final_table(:,2) = nanstd(table,1);
final_table(:,3) = min(table,[],1);
final_table(:,4) = max(table,[],1);


table_cols = names;

save('all_table.mat','table');

save('final_table.mat','final_table');
save('table_columns','table_cols');

    keyboard;


end


function params = local_func_find_params(filename)
    
load(filename);

[m,n] = size(lob);

S = datevec(myTime);  % change the format of time from datenum into date vector of matlal;

my_start_time_lob = find(S(:,4)==8,1,'first');  % we only start to count every day from 8:00 am; 
my_end_time_lob = find(S(:,4)<21,1,'last');

if isempty(my_start_time_lob)
    my_start_time_lob = 20000;
end

if isempty(my_start_time_lob)
    my_end_time_lob = m-20000;
end


my_start_time = find(arr(:,2)>= myTime(my_start_time_lob),1,'first');
my_end_time = find(arr(:,2) <= myTime(my_end_time_lob),1,'last');

% my_start_time = 20108;  % we ignore the very beginning period when the market almost has no trading. 
% % my_end_time =   2215165;  
%  my_end_time = 2215250; % we ignore the very last period when the market almost has no trading. 


start_time = datestr(arr(my_start_time,2),'dd-mm-yyyy HH:MM:SS.FFF'),
end_time = datestr(arr(my_end_time,2),'dd-mm-yyyy HH:MM:SS.FFF'),

T_trading = 600*(single(end_time(12))-single(start_time(12)))+60*(single(end_time(13))-single(start_time(13)))+ 10*(single(end_time(15))-single(start_time(15)))+ 1*(single(end_time(16))-single(start_time(16)))+(1/6)*(single(end_time(18))-single(start_time(18)))+(1/60)*(single(end_time(19))-single(start_time(19)));
T_trading = T_trading*600;
display(['time of trading = ', num2str(T_trading), ' 100 ms']);

Dmax = 10;

market_order_number = 0;
market_order_size = 0;

ask_market_order_number = 0;
ask_market_order_size = 0;

bid_market_order_number = 0;
bid_market_order_size = 0;

ask_limit_order_number = zeros(Dmax,1);
ask_limit_order_size = zeros(Dmax,1);

ask_cancel_order_number = zeros(Dmax,1);
ask_cancel_order_size = zeros(Dmax,1);

bid_limit_order_number = zeros(Dmax,1);
bid_limit_order_size = zeros(Dmax,1);

bid_cancel_order_number = zeros(Dmax,1);
bid_cancel_order_size = zeros(Dmax,1);


rare_cases = 0;    % the cases when trading decrease is not equal to the first change after trading;
rare_cases2 = 0;
% when there is a increase in volume, it comes from the tick which can
% easily be checked, as well as the number of orders;#
% when there is a decrease in volume, it might be two reason, a cancellation of
% outstanding order, or due to trade. If due to trade, 

k = 1;
j = 1;
n = 1;
market_order_index = ones(1,1);
theta_test = 0;
kkk=1;
jjj=1;

trade_count = 1;
trade_stack = nan(m,1);
limit_count = 1;
limit_stack = nan(m,1);
cancel_count = 1;
cancel_stack = nan(m,1);
for t = my_start_time: my_end_time
%     t = 840000+239;
    if deltaVol(t)>0
        pLevel = arr(t,1);
        nuOfOrders = arr(t,12);
        if(~(isnan(arr(t,3)) && isnan(arr(t,4))))
        % BID
        side = 2;
        price = arr(t,3);
        mySize = deltaVol(t);
        bid_limit_order_size(pLevel) = bid_limit_order_size(pLevel)+mySize;
        bid_limit_order_number(pLevel) = bid_limit_order_number(pLevel)+1;
        limit_stack(limit_count) = mySize;
        limit_count = limit_count+1;
        elseif(~(isnan(arr(t,5)) && isnan(arr(t,6))))
            % ASK
        side = 1;    
        price = arr(t,5);
        mySize = deltaVol(t); 
        ask_limit_order_size(pLevel) = ask_limit_order_size(pLevel)+mySize;
        ask_limit_order_number(pLevel) = ask_limit_order_number(pLevel)+1;
        limit_stack(limit_count) = mySize;
        limit_count = limit_count+1;        
        end
        
    elseif deltaVol(t)<0
        nuOfOrders = arr(t,12);
        pLevel = arr(t,1);
        mySize = abs(deltaVol(t));
        % there is two possibilities, due to trading or cancelling
        % if it is due to trading... we need a critertion;
        % SameTimestamp_arr = market_order_index(end)+find(((arr((market_order_index(end)+1):(t-1),2) == arr(t,2)).*(arr((market_order_index(end)+1):(t-1),1) == -1))); 
        
        
        if (arr(t,1) > 1)||((arr(t,1)==1)&&((~(isnan(arr(t,3)))&&(arr(t-1,11)==1))||((~(isnan(arr(t,5)))&&(arr(t-1,11)==2)))))
                % this is due to normal cancellations of order book;
                        if(~(isnan(arr(t,3)) && isnan(arr(t,4))))
                        % BID
                        side = 2;
                        price = arr(t,3);
                        mySize =  abs(deltaVol(t));
                        bid_cancel_order_size(pLevel) = bid_cancel_order_size(pLevel)+mySize;
                        bid_cancel_order_number(pLevel) = bid_cancel_order_number(pLevel)+1;
                        if pLevel ==1
                            theta_test(kkk) = mySize./(arr(t,4)+mySize);
                            theta_test2(kkk) = mySize;
                            theta_test3(kkk) = arr(t,4);
                            kkk=kkk+1;
                        end                       
                        cancel_stack(cancel_count) = mySize;
                        cancel_count = cancel_count+1;                           
                        
                        elseif(~(isnan(arr(t,5)) && isnan(arr(t,6))))
                        % ASK
                        side = 1;    
                        price = arr(t,5);
                        mySize =  abs(deltaVol(t));
                        ask_cancel_order_size(pLevel) = ask_cancel_order_size(pLevel)+mySize;
                        ask_cancel_order_number(pLevel) = ask_cancel_order_number(pLevel)+1;        
                        if pLevel ==1
                            theta_test(kkk) = mySize./(arr(t,6)+mySize);
                            theta_test2(kkk) = mySize;
                            theta_test3(kkk) = arr(t,6);
                            kkk=kkk+1;
                        end
                        cancel_stack(cancel_count) = mySize;
                        cancel_count = cancel_count+1;                            
                        end
        else
            
            trades_nearby_Endtime = (t-abs(deltaVol(t))-5000)-1+ find(isnan(deltaVol((t-abs(deltaVol(t))-5000):t-1)),1,'last');
            trades_nearby_Starttime = trades_nearby_Endtime; 
            i = 1;
            while i <= trades_nearby_Endtime - abs(deltaVol(t))
                if isnan(deltaVol(trades_nearby_Endtime-i))
                    i = i+1;
                else
                    trades_nearby_Starttime = trades_nearby_Endtime - i+1;
                    break;
                end
            end
            
            if (arr(t,1) ==1)&&(market_order_index(end)<trades_nearby_Endtime)&&(sum(arr(trades_nearby_Starttime:trades_nearby_Endtime,8))>=-deltaVol(t))     
            % determine if it is the nearest nagative volume, with tick =1, in the same side of trading and also explains the nearby trades volume;
                    market_order_number = market_order_number+ i;
                    market_order_size = -deltaVol(t) + market_order_size;
                    trade_stack(trade_count) = -deltaVol(t);
                    trade_count = trade_count + 1;
                    if arr(trades_nearby_Endtime,11)==1 % AggressorSide = buy ; otherwise sell side; buy side means ask side decrease
                        ask_market_order_number = ask_market_order_number + i;
                        ask_market_order_size = -deltaVol(t) + ask_market_order_size;
                    else
                        bid_market_order_number = bid_market_order_number + i;
                        bid_market_order_size = -deltaVol(t) + bid_market_order_size;
                    end
                    
                    market_order_index = [market_order_index,trades_nearby_Starttime:trades_nearby_Endtime];                
                    if sum(arr(trades_nearby_Starttime:trades_nearby_Endtime,8))~=-deltaVol(t)
                    	rare_cases = rare_cases+1;  % calculate the rare cases;
                        rare_case_index(k,1) = t;
                        rare_case_index(k,2) = -deltaVol(t) - sum(arr(trades_nearby_Starttime:trades_nearby_Endtime,8));
                        k = k+1;
                        mySize = abs(rare_case_index(k-1,2));
                        if rare_case_index(k-1,2)< 0  % there is a cancelling order occurring plus trading;
                            if(~(isnan(arr(t,3)) && isnan(arr(t,4))))
                            % BID
                            side = 2;
                            price = arr(t,3);
                            %mySize =  abs(deltaVol(t));
                            bid_cancel_order_size(pLevel) = bid_cancel_order_size(pLevel)+mySize;
                            bid_cancel_order_number(pLevel) = bid_cancel_order_number(pLevel)+1;
                            cancel_stack(cancel_count) = mySize;
                            cancel_count = cancel_count+1;                               
                            
                            if pLevel ==1
                                theta_test(kkk) = mySize./(arr(t,4)+mySize);
                                theta_test2(kkk) = mySize;
                                theta_test3(kkk) = arr(t,4);
                                kkk=kkk+1;
                            end
                            elseif(~(isnan(arr(t,5)) && isnan(arr(t,6))))
                            % ASK
                            side = 1;    
                            price = arr(t,5);
                            %mySize =  abs(deltaVol(t));
                            ask_cancel_order_size(pLevel) = ask_cancel_order_size(pLevel)+mySize;
                            ask_cancel_order_number(pLevel) = ask_cancel_order_number(pLevel)+1;
                            cancel_stack(cancel_count) = mySize;
                            cancel_count = cancel_count+1;                               
                            if pLevel ==1
                                theta_test(kkk) = mySize./(arr(t,6)+mySize);
                                theta_test2(kkk) = mySize;
                                theta_test3(kkk) = arr(t,6);
                                kkk=kkk+1;
                            end
                            end                            
                            
                        else    % there is a increasing order occurring plus trading;
                            if(~(isnan(arr(t,3)) && isnan(arr(t,4))))
                            % BID
                            side = 2;
                            price = arr(t,3);
                            %mySize =  abs(deltaVol(t));
                            limit_stack(limit_count) = mySize;
                            limit_count = limit_count+1;                                
                            bid_limit_order_size(pLevel) = bid_limit_order_size(pLevel)+mySize;
                            bid_limit_order_number(pLevel) = bid_limit_order_number(pLevel)+1;
                            elseif(~(isnan(arr(t,5)) && isnan(arr(t,6))))
                                % ASK
                            side = 1;    
                            price = arr(t,5);
                            %mySize =  abs(deltaVol(t));
                            ask_limit_order_size(pLevel) = ask_limit_order_size(pLevel)+mySize;
                            ask_limit_order_number(pLevel) = ask_limit_order_number(pLevel)+1;
                            limit_stack(limit_count) = mySize;
                            limit_count = limit_count+1;     
%                             ask_cancel_order_number(pLevel) = ask_cancel_order_number(pLevel)+1;        

                            end                        
                        end
                    end
            else  % normal cancellation orders;
                        if(~(isnan(arr(t,3)) && isnan(arr(t,4))))
                        % BID
                        side = 2;
                        price = arr(t,3);
                        mySize = abs(deltaVol(t));
                        bid_cancel_order_size(pLevel) = bid_cancel_order_size(pLevel)+mySize;
                        bid_cancel_order_number(pLevel) = bid_cancel_order_number(pLevel)+1;
                        cancel_stack(cancel_count) = mySize;
                        cancel_count = cancel_count+1;                           
                        if pLevel ==1
                            theta_test(kkk) = mySize./(arr(t,4)+mySize);
                            theta_test2(kkk) = mySize;
                            theta_test3(kkk) = arr(t,4);
                            kkk=kkk+1;
                        end
                        elseif(~(isnan(arr(t,5)) && isnan(arr(t,6))))
                        % ASK
                        side = 1;    
                        price = arr(t,5);
                        mySize =  abs(deltaVol(t));
                        ask_cancel_order_size(pLevel) = ask_cancel_order_size(pLevel)+mySize;
                        ask_cancel_order_number(pLevel) = ask_cancel_order_number(pLevel)+1; 
                        cancel_stack(cancel_count) = mySize;
                        cancel_count = cancel_count+1;                           
                        if pLevel ==1
                            theta_test(kkk) = mySize./(arr(t,6)+mySize);
                            theta_test2(kkk) = mySize;
                            theta_test3(kkk) = arr(t,6);
                            kkk=kkk+1;
                        end
                        end
            end
        end
    end
%     if mod(t,1000)==0
%     display([num2str(round((t-my_start_time)/(my_end_time-my_start_time)*100)),'%']),
%     end
end


%%

trades2 =find(arr(my_start_time:my_end_time,1) == -1);

market_order_number2 = length(trades2);
market_order_index = market_order_index(2:end);


limit_order_size = ask_limit_order_size + bid_limit_order_size;
limit_order_number = ask_limit_order_number + bid_limit_order_number;

cancel_order_size = ask_cancel_order_size + bid_cancel_order_size;
cancel_order_number = ask_cancel_order_number + bid_cancel_order_number;

% (cancel_order_size+ market_order_size)./limit_order_size,

avg_size_LO = sum(limit_order_size(1:9))./sum(limit_order_number(1:9)),
avg_size_LO(1) = 1;
avg_size_Cancel = sum(cancel_order_size(1:9))./sum(cancel_order_number(1:9)),
avg_size_Market = market_order_size./market_order_number,

avg_size_LO2 =  limit_order_size./limit_order_number;

% M = 200;
% distance = round((my_end_time - my_start_time)./M);


lob_start_time = my_start_time - length(find(arr(1:my_start_time-1,1)==-1));
lob_end_time = my_end_time - length(find(arr(1:my_end_time-1,1)==-1));


% size_lob = size(lob,1);
% size_arr = size(arr,1);
% 
% lob_start_time = round(my_start_time./size_arr.*size_lob);
% lob_end_time = round(my_end_time./size_arr.*size_lob);
% 
% distance = round((lob_end_time - lob_start_time)./M);
 estimate_Q_range = lob_start_time:1:lob_end_time;


% catlob=cat(4,lob{estimate_Q_range});
% 
% Q_a = squeeze(sum(catlob(1,2,:,:),4)./avg_size_LO./length(estimate_Q_range));
% Q_b = squeeze(sum(catlob(2,2,:,:),4)./avg_size_LO./length(estimate_Q_range));
% 
% Q = (Q_a+Q_b)./2;
% bar(avg_size_LO);
% title('average
%%
% [a1, b1] = cce_LUCA_powerLawFit(lamda_est(1:5), [1:5]');
% 
% figure();
% % title('Power law fit of Limit orders rate');
% 
% 
% % bar(limit_order_number./T_trading);
% % hold on;bar(a1.*[1:5].^b1,'r');
% lamda_fit = [limit_order_number./T_trading,a1.*[1:10]'.^b1];
% 
% lamda_fit(6:end,2) = nan;
% bar(lamda_fit);
% xlabel('distance from opposite quote');
% ylabel('\lambda');
% legend('DATA','MODEL');
% % export_fig 'fig_powerlawfit_lambda.pdf' -pdf -transparent;
% 
% 
% figure();
% bar(cancel_order_number./T_trading./Q.*avg_size_Cancel./avg_size_LO);
% title('Cancellation rates');
% xlabel('distance from opposite quote');


%%

% lamda_est = limit_order_size./T_trading./avg_size_LO;
% seita_est =
% cancel_order_number./T_trading./Q.*avg_size_Cancel./avg_size_LO(1);  %
% Cont 2010 approach
% u = market_order_number./T_trading*avg_size_Market./avg_size_LO2(1);
% u = market_order_number./T_trading*avg_size_Market./avg_size_LO;

% Our own approach (similar to Cont 2013)
% seita_est = cancel_order_number./T_trading.*avg_size_Cancel./avg_size_LO;

% save('parameters_estimation_2013Approach','lamda_est','seita_est','u','Q','T_trading','limit_order_number','limit_order_size','cancel_order_number','theta_test','cancel_order_size','avg_size_Cancel','avg_size_LO','avg_size_Market','market_order_number','market_order_size');

% (u+seita_est(1))/lamda_est(1),
% 
% [a1, b1] = cce_LUCA_powerLawFit(lamda_est(1:5), [1:5]');
% [a2, b2] = cce_LUCA_powerLawFit(seita_est(1:5), [1:5]');

% plot([1:5],a2.*[1:5].^b2); hold on; plot([1:5],seita_est(1:5),'*r');

% GGG = [theta_test;theta_test2;theta_test3];
% 
% sort_cancel = sortrows(GGG',2)';
% 
% const_index1 = find(sort_cancel(2,:)<=5);
% seita_const1 = sum(sort_cancel(2,const_index1))./market_order_size*u;
% 
% sort_cancel_left = sortrows(sort_cancel(:,const_index1+1:end)',1)';
% prop_index = [find(log10(sort(sort_cancel(1,const_index1+1:end)))>=-2.5,1,'first'):find(log10(sort(sort_cancel(1,const_index1+1:end)))<=-1,1,'last')];
% 
% prop_const = mean(sort_cancel_left(1,prop_index));
% 
% seita_const2 = sum(sort_cancel_left(2,[1:prop_index(1)-1,prop_index(end)+1:end]))./market_order_size*u;
% 
% seita_const_all = sum(sort_cancel(2,:))./market_order_size*u;
% 
% sort_cancel_byProp = sortrows(GGG',1)';
% 
% seita_const3 = sum(sort_cancel_byProp(2,1:3.594*10^5))./market_order_size*u;
% 
% %option 1
% u+seita_const2+seita_const1,lamda_est(1),
% 
% %option 2,
% seita_const_all+u,lamda_est(1),

%option 3,

% theta_all2 = struct('lambda_a',lambda_a,'mu_a',mu_a,'theta_a',theta_a,'lambda_b',lambda_b,'mu_b',mu_b,'theta_b',theta_b);
% lambda_a = ask_limit_order_size(1)./T_trading./avg_size_LO2(1);
% lambda_b = bid_limit_order_size(1)./T_trading./avg_size_LO2(1);
% mu_a = ask_market_order_size./T_trading./avg_size_LO2(1);  % suppose balanced here; needs to be recalculated later;
% mu_b = bid_market_order_size./T_trading./avg_size_LO2(1);
% theta_a = ask_cancel_order_size(1)./T_trading./avg_size_LO2(1);
% theta_b = bid_cancel_order_size(1)./T_trading./avg_size_LO2(1);

% plot(log10(sort_cancel(2,3.1415*10^5:end)));
% figure();
% plot(log10(sort_cancel(1,3.1415*10^5:end)));
% figure();
% plot(log10((sort(sort_cancel(1,3.1415*10^5:end)))));
 %save('parameters_estimation','lamda_est','seita_est','u','Q','T_trading','limit_order_number','limit_order_size','cancel_order_number','theta_test','cancel_order_size','avg_size_Cancel','avg_size_LO','avg_size_Market','market_order_number','market_order_size');

%%
% to produce another set of parameters for un-balanced order flow for Cont
% 2013 model.


% lamda_est = limit_order_size./T_trading./avg_size_LO;
% u = ask_market_order_size/T_trading*avg_size_Market./avg_size_LO;
% seita_est = ask_cancel_order_number./T_trading.*avg_size_Cancel./avg_size_LO;


% ask_dec_size = ask_market_order_size+ask_cancel_order_size;
% ask_inc_size = ask_limit_order_size;
% bid_dec_size = bid_market_order_size + bid_cancel_order_size;
% bid_inc_size = bid_limit_order_size;

% p_20 = ask_dec_size./T_trading./avg_size_LO;
% p_10 = ask_inc_size./T_trading./avg_size_LO;
% p_01 = bid_inc_size./T_trading./avg_size_LO;
% p_02 = bid_dec_size./T_trading./avg_size_LO;
% save('parameters_estimation_2013Cont_FTApproach','p_01','p_02','p_10','p_20');



params.totalVolume = bid_market_order_size + ask_market_order_size;
params.tradeSize_mean  = params.totalVolume./(ask_market_order_number+bid_market_order_number);
params.tradeSize_std  = nanstd(trade_stack);
params.tradeSize_max = max(trade_stack);
params.tradeSize_min = min(trade_stack);

params.numOrderTop_mean = nanmean(lob(my_start_time_lob:my_end_time_lob,5),1) + nanmean(lob(my_start_time_lob:my_end_time_lob,6),1);
params.numOrderTop_std = nanstd(nansum(lob(my_start_time_lob:my_end_time_lob,5:6),2));
params.numOrderTop_max = max(nansum(lob(my_start_time_lob:my_end_time_lob,5:6),2));
params.numOrderTop_min = min(nansum(lob(my_start_time_lob:my_end_time_lob,5:6),2));

% params.numOrderTop_max = max(max(lob(my_start_time_lob:my_end_time_lob,5),[],1),max(lob(my_start_time_lob:my_end_time_lob,6),[],1));
params.numOrderDeep_mean = nansum(nanmean(lob(my_start_time_lob:my_end_time_lob,[11:6:53]),1) + nanmean(lob(my_start_time_lob:my_end_time_lob,[12:6:54]),1));
params.numOrderDeep_std = nanstd(nansum(lob(my_start_time_lob:my_end_time_lob,[11:6:53,12:6:54]),2),[],1);
params.numOrderDeep_max = max(nansum(lob(my_start_time_lob:my_end_time_lob,[11:6:53,12:6:54]),2),[],1);
params.numOrderDeep_min = min(nansum(lob(my_start_time_lob:my_end_time_lob,[11:6:53,12:6:54]),2),[],1);

marketMidPrice = nansum(lob(my_start_time_lob:my_end_time_lob,1:2),2)./2;
marketReturns = diff(marketMidPrice);
marketReturns = [0; marketReturns];

market_bidask_spread = abs(lob(my_start_time_lob:my_end_time_lob,1) - lob(my_start_time_lob:my_end_time_lob,2));

myTime_selected = myTime(my_start_time_lob:my_end_time_lob);

params.numPriceChange = sum(marketReturns~=0);

time_when_market_price_change = myTime_selected(marketReturns~=0);

time_interval = diff(time_when_market_price_change).*3600*24; % calculate the time interval length between market price changes; units in seconds;


params.avgPriceMoveTime = nanmean(time_interval);  % this makes the units to be seconds. 
params.stdPriceMoveTime = nanstd(time_interval);
params.maxPriceMoveTime = max(time_interval);
params.minPriceMoveTime = min(time_interval);

params.meanLimitOrderSize = nanmean(limit_stack);
params.stdLimitOrderSize = nanstd(limit_stack);
params.minLimitOrderSize = min(limit_stack);
params.maxLimitOrderSize = max(limit_stack);

params.meanCancelOrderSize = nanmean(cancel_stack);
params.stdCancelOrderSize = nanstd(cancel_stack);
params.minCancelOrderSize = min(cancel_stack);
params.maxCancelOrderSize = max(cancel_stack);


params.meanBidAskSpread = nanmean(market_bidask_spread);
params.stdBidAskSpread = nanstd(market_bidask_spread);
params.minBidAskSpread = min(market_bidask_spread);
params.maxBidAskSpread = max(market_bidask_spread);

ts = 10;  % units in seconds. 
[params2,myTimeVector2] = cce_GZ2014_find_param_alldepth(arr,deltaVol);

state_display = ['parameter found'],
% ts = 10;
% myTimeVec = [myTimeVector2(1),myTimeVector2(end)];
[time_syn,~] = fun_grid_time([myTimeVector2(1),myTimeVector2(end)],ts);          % this generates a time grid with ts as time step starting from 22:45 to next day 21:15;
params_time = function_syc(params2,myTimeVector2,time_syn');
params_time_diff = [params_time(1,:);diff(params_time)];
params_time_diff(isnan(params_time_diff)) =0;
params3 = params_time_diff;   % this is the resampled time vector; 
% [params,myTimeVec,lob] = cce_GZ2014a_simulation_fast_online_find_params_resample(filename,ts);

hit_level =9;
order_flow(:,1) = nansum(params3(:,1:1+hit_level-1),2); % ASK flow increase in sum; 
order_flow(:,2) = abs(nansum(params3(:,11:(11+hit_level-1)),2) - params3(:,41)); % ASK flow cancelled decrease in sum; 
order_flow(:,3) = nansum(params3(:,21:(21+hit_level-1)),2); % BID flow increase in sum; 
order_flow(:,4) = abs(nansum(params3(:,31:(31+hit_level-1)),2) - params3(:,42)); % BID flow cancelled decrease in sum; 
order_flow(:,5) = params3(:,41);  % Market buy order;
order_flow(:,6) = params3(:,42);  % Market sell order; 

params.meanLimitOrderRate = nanmean(order_flow(:,1) + order_flow(:,3))./10;
params.stdLimitOrderRate = nanstd(order_flow(:,1)./10 + order_flow(:,3)./10);
params.minLimitOrderRate = min(order_flow(:,1)./10 + order_flow(:,3)./10);
params.maxLimitOrderRate = max(order_flow(:,1)./10 + order_flow(:,3)./10);

params.meanCancelOrderRate = nanmean(order_flow(:,2) + order_flow(:,4))./10;
params.stdCancelOrderRate = nanstd(order_flow(:,2)./10 + order_flow(:,4)./10);
params.minCancelOrderRate = min(order_flow(:,2)./10 + order_flow(:,4)./10);
params.maxCancelOrderRate = max(order_flow(:,2)./10 + order_flow(:,4)./10);


params.meanMarketOrder = nanmean(order_flow(:,5)./10+ order_flow(:,6)./10);
params.stdMarketOrder = nanstd(order_flow(:,5)./10+ order_flow(:,6)./10);
params.minMarketOrder = min(order_flow(:,5)./10+ order_flow(:,6)./10);
params.maxMarketOrder = max(order_flow(:,5)./10+ order_flow(:,6)./10);


end