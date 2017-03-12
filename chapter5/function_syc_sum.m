function ts3Data = function_syc_sum(data1,time1,time2)

% [data1 time1] = cce_split(ts1);
% 
% [data2 time2] = cce_split(ts2);
% 
% clear data2;
% data1 = bid_dec_cum;

nLen  = length(time2);
nLen2 = length(time1);
ts3Data = NaN(nLen,size(data1,2));

% idxTemp = -999;
idxTemp2 = 1;
idxshift = 20000;
idx2 = 1;
%%

% if(strcmpi(type, 'default'))

    %the standard one we use when fitting mkt returns

   

    for i = 1: nLen
        
        if idxTemp2+idxshift > nLen2
            idxshift = nLen2 - idxTemp2;
        end

        idx = find(time1(idxTemp2:(idxTemp2+idxshift)) <= time2(i), 1,'last');
        
        if ~isempty(idx)
            idx = idx + idxTemp2 - 1;
        end
        
        if isempty(idx)
            idx = find(time1 <= time2(i), 1,'last');
        end
        
        if(isempty(idx))

            ts3Data(i,:) = 0;

        else
            if idx2~=idx
                ts3Data(i,:) = sum(data1((idx2+1):idx,:));
            else
                ts3Data(i,:) = 0;
            end

        end       

        %be careful not to just propagate values through. eg assign returns to

        %a timestamp where they did not happen etc

%         if(~isequal(idx, idxTemp))

            %we have found a new value for idx
            idx2 = idx;
            idxTemp = idx;
            idxshift = min(max(idxshift, (idxTemp - idxTemp2)*2),100000);
            if isempty(idxshift) 
                idxshift = 20000;
            end
            idxTemp2 = idx;

%         else
% 
%             %ensure that NaN is inserted into the output sequence
%             idx = [];
% 
%             idx = [];
% 
%         end



        clear idx;
        
        if mod(i,10000) == 0
            (i./nLen)*100,
        end

    end

end


%    function [ret,current_flag] = sum_return(last_flag)
%        
%         
%     % Each time this function searches after last index of market return timeseries to find the next index of timestamps of market
%     % returns which is larger than the current timestamp (DATA_TIME_lag(n), changes as the main program runs) of trading signal; then
% 
%         current_flag = find(DATA_TIME(last_flag:min(last_flag+20000,length(DATA_TIME)))> DATA_TIME_lag(n),1,'first');
%         
%         if isempty(current_flag)
%               current_flag = find(DATA_TIME(last_flag:end)> DATA_TIME_lag(n),1,'first');
%         end
%         current_flag = current_flag + last_flag -2;
%         
%         if current_flag == last_flag
%             ret = 0;
%         else
%             if last_flag == 1;
%                 last_flag = 0;
%             end
%             ret = sum(Market_return((last_flag+1):current_flag));      
%             % all the market returns happened between these two indices are the
%             % returns we have got between the trading signals. Stored in
%             % ret.
%         end
%     end
