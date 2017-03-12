function [params,myTimeVec] = cce_GZ2014_find_param_alldepth(arr,deltaVol)

% find a start time to process; and an end time to end;
% default start time  =  22:45 previous day; 
% default end time = 21:15 next day;
% if data date is shorter than this interval, then process the nearest
% hours, e.g. 22:00 to 21:00, etc...

start_date = datevec(arr(1,2));
end_date = datevec(arr(end,2));


if (start_date(4)*60+start_date(5)<22*60+45)&&(start_date(3) ~= end_date(3))
    start_date_vec = start_date;
    start_date_vec(4) = 22;
    start_date_vec(5) = 45;
    start_date_vec(6) = 0;
    start_date_num = datenum(start_date_vec);
else
    start_date_num = ceil(arr(1,2)*24)/24;
end

if (end_date(4)*60+end_date(5)<22*60+45)
    end_date_vec = end_date;
    end_date_vec(4) = 21;
    end_date_vec(5) = 15;
    end_date_vec(6) = 20;
    end_date_num = datenum(end_date_vec);
else
    end_date_num = fix(arr(end,2)*24)/24;
end

my_start_time = find(arr(:,2)>start_date_num,1,'first');
my_end_time = find(arr(:,2)<end_date_num,1,'last');

%%
% start_time = datestr(arr(my_start_time,2),'dd-mm-yyyy HH:MM:SS.FFF'),
% end_time = datestr(arr(my_end_time,2),'dd-mm-yyyy HH:MM:SS.FFF'),

T_trading = (arr(my_end_time,2) - arr(my_start_time,2))*24*36000;
% T_trading_time = 600*(single(end_time(12))-single(start_time(12)))+60*(single(end_time(13))-single(start_time(13)))+ 10*(single(end_time(15))-single(start_time(15)))+ 1*(single(end_time(16))-single(start_time(16)))+(1/6)*(single(end_time(18))-single(start_time(18)))+(1/60)*(single(end_time(19))-single(start_time(19)));
% T_trading = T_trading_time*600;   % set t_trading into unit of 100ms;
% display(['time of trading = ', num2str(T_trading), ' 100 ms']);
myTimeVec = arr(my_start_time:my_end_time,2);
NN = length(myTimeVec);
ask_inc_order_size = nan(NN,10);
bid_inc_order_size = nan(NN,10);
ask_dec_order_size = nan(NN,10);
bid_dec_order_size = nan(NN,10);
bid_market_order_size_all = nan(NN,1);
ask_market_order_size_all = nan(NN,1);
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
for t = my_start_time: my_end_time
% for t = my_start_time: (my_start_time+50000)

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
        elseif(~(isnan(arr(t,5)) && isnan(arr(t,6))))
            % ASK
        side = 1;    
        price = arr(t,5);
        mySize = deltaVol(t); 
        ask_limit_order_size(pLevel) = ask_limit_order_size(pLevel)+mySize;
        ask_limit_order_number(pLevel) = ask_limit_order_number(pLevel)+1;        
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
%                         if pLevel ==1
%                             theta_test(kkk) = mySize./(arr(t,4)+mySize);
%                             theta_test2(kkk) = mySize;
%                             theta_test3(kkk) = arr(t,4);
%                             kkk=kkk+1;
%                         end                       
                        
                        
                        elseif(~(isnan(arr(t,5)) && isnan(arr(t,6))))
                        % ASK
                        side = 1;    
                        price = arr(t,5);
                        mySize =  abs(deltaVol(t));
                        ask_cancel_order_size(pLevel) = ask_cancel_order_size(pLevel)+mySize;
                        ask_cancel_order_number(pLevel) = ask_cancel_order_number(pLevel)+1;        
%                         if pLevel ==1
%                             theta_test(kkk) = mySize./(arr(t,6)+mySize);
%                             theta_test2(kkk) = mySize;
%                             theta_test3(kkk) = arr(t,6);
%                             kkk=kkk+1;
%                         end
                        
                        end
        else
            
            
            if (t-abs(deltaVol(t))-5000)<1
                continue;
            else
            
                trades_nearby_Endtime = (t-abs(deltaVol(t))-5000)-1+ find(isnan(deltaVol((t-abs(deltaVol(t))-5000):t-1)),1,'last');
                if isempty(trades_nearby_Endtime)
                    trades_nearby_Endtime = find(isnan(deltaVol(1:t-1)),1,'last');
                end
            end
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
%                             if pLevel ==1
%                                 theta_test(kkk) = mySize./(arr(t,4)+mySize);
%                                 theta_test2(kkk) = mySize;
%                                 theta_test3(kkk) = arr(t,4);
%                                 kkk=kkk+1;
%                             end
                            elseif(~(isnan(arr(t,5)) && isnan(arr(t,6))))
                            % ASK
                            side = 1;    
                            price = arr(t,5);
                            %mySize =  abs(deltaVol(t));
                            ask_cancel_order_size(pLevel) = ask_cancel_order_size(pLevel)+mySize;
                            ask_cancel_order_number(pLevel) = ask_cancel_order_number(pLevel)+1;        
%                             if pLevel ==1
%                                 theta_test(kkk) = mySize./(arr(t,6)+mySize);
%                                 theta_test2(kkk) = mySize;
%                                 theta_test3(kkk) = arr(t,6);
%                                 kkk=kkk+1;
%                             end
                            end                            
                            
                        else    % there is a increasing order occurring plus trading;
                            if(~(isnan(arr(t,3)) && isnan(arr(t,4))))
                            % BID
                            side = 2;
                            price = arr(t,3);
                            %mySize =  abs(deltaVol(t));
                            bid_limit_order_size(pLevel) = bid_limit_order_size(pLevel)+mySize;
                            bid_limit_order_number(pLevel) = bid_limit_order_number(pLevel)+1;
                            elseif(~(isnan(arr(t,5)) && isnan(arr(t,6))))
                                % ASK
                            side = 1;    
                            price = arr(t,5);
                            %mySize =  abs(deltaVol(t));
                            ask_limit_order_size(pLevel) = ask_limit_order_size(pLevel)+mySize;
                            ask_limit_order_number(pLevel) = ask_limit_order_number(pLevel)+1;        
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
%                         if pLevel ==1
%                             theta_test(kkk) = mySize./(arr(t,4)+mySize);
%                             theta_test2(kkk) = mySize;
%                             theta_test3(kkk) = arr(t,4);
%                             kkk=kkk+1;
%                         end
                        elseif(~(isnan(arr(t,5)) && isnan(arr(t,6))))
                        % ASK
                        side = 1;    
                        price = arr(t,5);
                        mySize =  abs(deltaVol(t));
                        ask_cancel_order_size(pLevel) = ask_cancel_order_size(pLevel)+mySize;
                        ask_cancel_order_number(pLevel) = ask_cancel_order_number(pLevel)+1;        
%                         if pLevel ==1
%                             theta_test(kkk) = mySize./(arr(t,6)+mySize);
%                             theta_test2(kkk) = mySize;
%                             theta_test3(kkk) = arr(t,6);
%                             kkk=kkk+1;
%                         end
                        end
            end
            
%             end
        end
    end
    ask_dec_order_size(t-my_start_time+1,1)=  ask_market_order_size(1)+ask_cancel_order_size(1);
    bid_dec_order_size(t-my_start_time+1,1)=  bid_market_order_size(1)+bid_cancel_order_size(1);
%   bid_market_order_size = nan(NN,1);
%   ask_market_order_size = nan(NN,1);
    ask_market_order_size_all(t-my_start_time+1,1) = ask_market_order_size(1);
    bid_market_order_size_all(t-my_start_time+1,1) = bid_market_order_size(1);
    ask_inc_order_size(t-my_start_time+1,1)=  ask_limit_order_size(1);
    bid_inc_order_size(t-my_start_time+1,1)=  bid_limit_order_size(1);
    
% here we add to include all other depth-data.    
    
%     ask_dec_order_size(t-my_start_time+1,2:end)=  ask_cancel_order_size(2:end);
%     bid_dec_order_size(t-my_start_time+1,2:end)=  bid_cancel_order_size(2:end);
%     ask_inc_order_size(t-my_start_time+1,2:end)=  ask_limit_order_size(2:end);
%     bid_inc_order_size(t-my_start_time+1,2:end)=  bid_limit_order_size(2:end);
    ask_dec_order_size(t-my_start_time+1,2:end-1)=  ask_cancel_order_size(2:end-1);
    bid_dec_order_size(t-my_start_time+1,2:end-1)=  bid_cancel_order_size(2:end-1);
    ask_inc_order_size(t-my_start_time+1,2:end-1)=  ask_limit_order_size(2:end-1);
    bid_inc_order_size(t-my_start_time+1,2:end-1)=  bid_limit_order_size(2:end-1);   
%     if mod(t,1000)==0
%     display([num2str(round((t-my_start_time)/(my_end_time-my_start_time)*100)),'%']),
%     end
end

params = [ask_inc_order_size,ask_dec_order_size,bid_inc_order_size,bid_dec_order_size,ask_market_order_size_all,bid_market_order_size_all];


% ask_dec_size = ask_market_order_size(1)+ask_cancel_order_size(1);
% ask_inc_size = ask_limit_order_size(1);
% bid_dec_size = bid_market_order_size(1) + bid_cancel_order_size(1);
% bid_inc_size = bid_limit_order_size(1);
% 
% p20 = ask_dec_size./T_trading;
% p10 = ask_inc_size./T_trading;
% p01 = bid_inc_size./T_trading;
% p02 = bid_dec_size./T_trading;
% 
% params = [p10,p20,p01,p02];   % ask side inc, ask side dec, bid side inc, bid side dec;

end
