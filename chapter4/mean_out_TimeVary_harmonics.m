function [T_out2,std_T_out2] = mean_out_TimeVary_harmonics(w_new,x_p,data_pointer);   


        [m1,m2,m3] = size(x_p);
        ss =0;
        ss2 = zeros(m1,1);
        std_T_out2 = 0;
        for k=1:m1
                N_total_period = length(x_p(k,:,end));
                x_p_cumPeriod = zeros(N_total_period,1);
%                 for i = 1:N_total_period
%                     x_p_cumPeriod(i) = sum(x_p(k,1:i,end));
%                 end
                x_p_cumPeriod = cumsum(x_p(k,:,end));
                P_left = max(find(x_p_cumPeriod - min(data_pointer)<0));  % find the i_th period in the left of data_pointer;
                P_right = min(find(x_p_cumPeriod - max(data_pointer)>0)); % find the j_th period in the right of data_pointer;
            if (isempty(P_left)==1)||(P_right<=1)
                x_p_out = mean(x_p,1);        
                T_out2 = x_p_out(1,1,end);
                std_T_out2 = 0;
                return;
            else
                ss = w_new(k).*x_p(k,P_right,end)+ss;
                ss2(k) = x_p(k,P_right,end);
            end
        end
        
        T_out2 = ss;
        std_T_out2 = sqrt(dot((ss2-ss).^2,w_new));
end