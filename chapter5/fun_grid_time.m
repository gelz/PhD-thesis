function [time_syn,grid_index] =fun_grid_time(time_vector,ts)

% A function that generate another time vector with uniform sampling time
% interevals specified on ts (in seconds),  start at the minimum and
% maximum of input time_vector, with other time located around nearest
% integrated hours.


ts1 = find_start_hours(time_vector);  % unit in day for time_vector;
ts2 = find_end_hours(time_vector);

ts_step = ts./24/3600;

time_syn = [time_vector(1),fliplr([ts1-ts_step:-ts_step:time_vector(1)]),ts1:ts_step:ts2,(ts2+ts_step):ts_step:time_vector(end),time_vector(end)];
grid_index = [2+length([ts1-ts_step:-ts_step:time_vector(1)]),find(time_syn==ts2)];
end

function ts1 = find_start_hours(time_vector)

dVec = datevec(time_vector(1));
% tindex2 = find((dVec(:,3)-3).*24+dVec(:,4)+dVec(:,5)./60+dVec(:,6)./3600>=1,1,'first');  % find out the 24:00 date index;

dVec2 = dVec;
dVec2(:,4) = dVec2(:,4)+1;
dVec2(:,5) = 0;
dVec2(:,6) = 0;
if dVec2(:,4) >23
    dVec2(:,4) = dVec2(:,4)-24;
    dVec2(:,3) = dVec2(:,3)+1;
end
ts1 = datenum(dVec2);
end

function ts2 = find_end_hours(time_vector)

dVec = datevec(time_vector(end));
% tindex2 = find((dVec(:,3)-3).*24+dVec(:,4)+dVec(:,5)./60+dVec(:,6)./3600>=1,1,'first');  % find out the 24:00 date index;
dVec2 = dVec;
dVec2(:,5) = 0;
dVec2(:,6) = 0;
ts2 = datenum(dVec2);

end