function  model_data = model_data_ARTimeVaryingFourierSeries_sin(x_p,data_pointer,data_sample,B_size);

    K_order = 40;
    p_order = 11;


    if min(data_pointer) ==1

        N_length = max(data_pointer);
        T0 = round(x_p(1,1,end));
        unit_freq = 2*pi/T0;
        unit_phase = unit_freq *(1:N_length);

    for k = 1 : (K_order/2)
        c(:,k) = cos((k-1)*unit_phase);
    end
    for k = ((K_order/2)+1) : K_order
        c(:,k) = sin((k-K_order/2)*unit_phase);
    end
        x_p_est = squeeze(x_p(1,1,1:end-1));
        harmonics_sum = sum(c*x_p_est(p_order+1:end),2);
        a= [1,x_p_est(1:p_order)'];
        output = filter(1,a,harmonics_sum);
        model_data = output(data_pointer)';
        return;
    end
    
    N_total_period = length(x_p(1,:,1));
%     x_p_cumPeriod = zeros(N_total_period,1);
%     for i = 1:N_total_period
%         x_p_cumPeriod(i) = sum(x_p(1,1:i,end));
%     end
    x_p_cumPeriod = cumsum((x_p(1,:,end)));
     
    P_left = max(find(x_p_cumPeriod - min(data_pointer)<0));  % find the i_th period in the left of data_pointer;
    P_right = min(find(x_p_cumPeriod - max(data_pointer)>0)); % find the j_th period in the right of data_pointer;
    if isempty(P_right)==1
        P_right = N_total_period;
    end
    if isempty(P_left)==0
        P_length = P_right - P_left;
    else
        N_length = max(data_pointer);
        T0 = round(x_p(1,1,end));
        unit_freq = 2*pi/T0;
        unit_phase = unit_freq *(1:N_length);

%     for k = 1 : (K_order/2)
% %         c(:,k) = cos((k-1)*unit_phase);
%     end
    c(:,1:(K_order/2)) = cos(unit_phase'*(0:(K_order/2-1)));
    c(:,(((K_order/2)+1) : K_order)) = sin(unit_phase'*(1:(K_order/2)));
%     for k = ((K_order/2)+1) : K_order
%         c(:,k) = sin((k-K_order/2)*unit_phase);
%     end
        x_p_est = reshape((x_p(1,1,1:end-1)),length(x_p(1,1,1:end-1)),1);
        harmonics_sum = sum(c*x_p_est(p_order+1:end),2);
        a= [1,x_p_est(1:p_order)'];
        output = filter(1,a,harmonics_sum);
        model_data = output(data_pointer)';
        return;
    end
    
    for  i = P_left:P_right-1


        T0 = round(x_p(1,i+1,end));
        unit_freq = 2*pi/T0;
        N_length = T0;
        unit_phase = unit_freq *(1:N_length);
        c = zeros(T0,K_order);
        c(:,1:(K_order/2)) = cos(unit_phase'*(0:(K_order/2-1)));
        c(:,(((K_order/2)+1) : K_order)) = sin(unit_phase'*(1:(K_order/2)));
%     for k = 1 : (K_order/2)
%         c(:,k) = cos((k-1)*unit_phase);
%     end
%     for k = ((K_order/2)+1) : K_order
%         c(:,k) = sin((k-K_order/2)*unit_phase);
%     end
       x_p_est = reshape((x_p(1,i+1,1:end-1)),length(x_p(1,i+1,1:end-1)),1);

%         x_p_est = (x_p(1,i+1,1:end-1));
        harmonics_sum = sum(c*x_p_est(p_order+1:end),2);
        a= [1,x_p_est(1:p_order)'];
        y_model(x_p_cumPeriod(i)+1:x_p_cumPeriod(i+1)) = filter(1,a,harmonics_sum);

    end
    
    model_data = y_model(data_pointer);


end


% function y = sum_harmonics(x_p(k,:,:))
%     
% 
% 
% 
% end




