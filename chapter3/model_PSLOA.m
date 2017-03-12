function [model_data,D_store] = model_PSLOA(x_p,data_pointer,data_sample,B_size,s0,D_store,k,seita_0)

%     model_data =  model_PSLOA(x_p(k,:,:),data_pointer,data_sample)

% for t = 1 : mod(length(data_sample),B_size)
%     data_pointer = ((t-1)*B_size+1):t*B_size;

    if min(data_pointer) ==1
        
        model_data = s0(data_pointer)';
        D_store(k,data_pointer) = s0(data_pointer);
        return;
    end
    
    N_total_period = length(x_p(1,:,1));
    x_p_cumPeriod = zeros(N_total_period,1);
    for i = 1:N_total_period
        x_p_cumPeriod(i) = sum(x_p(1,1:i,1));
    end
     
    P_left = max(find(x_p_cumPeriod - min(data_pointer)<0));  % find the i_th period in the left of data_pointer;
    P_right = min(find(x_p_cumPeriod - max(data_pointer)>0)); % find the j_th period in the right of data_pointer;
    if isempty(P_right)==1
        P_right = N_total_period;
    end
    if isempty(P_left)==0
        P_length = P_right - P_left;
    else
        model_data = s0(data_pointer)';
        D_store(k,data_pointer) = model_data;
        return;
    end
    
    for  i = P_left:P_right-1
        if i == P_left
            if P_left ==1
                s = x_p(1,(i+1),2)*D_store(k,1:x_p_cumPeriod(i));
            else
                s = x_p(1,(i+1),2)*D_store(k,x_p_cumPeriod(i-1)+1:x_p_cumPeriod(i));
            end
            if x_p(1,i+1,1)>=length(s)
                y_model(x_p_cumPeriod(i)+1:x_p_cumPeriod(i+1)) = [s,zeros(1,x_p(1,i+1,1)-length(s))];
            else
                y_model(x_p_cumPeriod(i)+1:x_p_cumPeriod(i+1)) = s(1:(x_p_cumPeriod(i+1)-x_p_cumPeriod(i)));
            end
        else
            s = x_p(1,(i+1),2)*y_model(x_p_cumPeriod(i-1)+1:x_p_cumPeriod(i));
            if x_p(1,i+1,1)>=length(s)
                y_model(x_p_cumPeriod(i)+1:x_p_cumPeriod(i+1)) = [s,zeros(1,x_p(1,i+1,1)-length(s))];
            else
                y_model(x_p_cumPeriod(i)+1:x_p_cumPeriod(i+1)) = s(1:(x_p_cumPeriod(i+1)-x_p_cumPeriod(i)));
            end
        end
    end
    [M1,N1] = size(y_model(x_p_cumPeriod(i)+1:x_p_cumPeriod(i+1)));
    y_model(x_p_cumPeriod(i)+1:x_p_cumPeriod(i+1)) = y_model(x_p_cumPeriod(i)+1:x_p_cumPeriod(i+1)) + seita_0(3)*randn(M1,N1);
    model_data = y_model(data_pointer);
    D_store(k,x_p_cumPeriod(i)+1:x_p_cumPeriod(i+1)) = y_model(x_p_cumPeriod(i)+1:x_p_cumPeriod(i+1));

end