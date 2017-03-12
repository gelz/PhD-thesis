function ts3 = fit_fast(ts1,ts2, type)

if(nargin<3)

    type = 'default';

end

if(isempty(ts1))

    ts3 = ts1;

    return;

end

%% Are they already fitted?

if(isequal(ts1.time, ts2.time))

    ts3 = ts1;

    return;

end

%%

 

%%

[data1 time1] = cce_split(ts1);

[data2 time2] = cce_split(ts2);

clear data2;

nLen  = length(time2);
nLen2 = length(time1);
ts3Data = NaN(nLen,size(data1,2));

idxTemp = -999;
idxTemp2 = 1;
idxshift = 20000;
%%

if(strcmpi(type, 'default'))

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
        
        

        %be careful not to just propagate values through. eg assign returns to

        %a timestamp where they did not happen etc

%         if(~isequal(idx, idxTemp))

            %we have found a new value for idx

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

        if(isempty(idx))

            ts3Data(i,:) = NaN;

        else

            ts3Data(i,:) = data1(idx,:);

        end

        clear idx;
        
%         if mod(i,10000) == 0
%             (i./nLen)*100,
%         end

    end

   

else

    dbstop if error;

    error('type not recognized')

end

 

%%

ts3 = timeseries(ts3Data, time2, 'Name', ts1.name); 

 

end

 