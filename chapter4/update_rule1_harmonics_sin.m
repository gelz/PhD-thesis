function  y = update_rule1_harmonics_sin(x_p,data_pointer,seita_0)

%         x_p(k,:,:) = update_rule1_harmonics(x_p(k,:,:),data_pointer,seita_0);
% var_period = 2;
% var_amplitude = 0.05;
% var_amplitude_ar = 0.05;
% 
% hyper_par = [var_period, var_amplitude,var_amplitude_ar];
p_order = 11;
K_order = 40;
T_var = seita_0(1);

        [m1,m2,m3] = size(x_p);
        y=x_p;
        for k=1:m1
                N_total_period = length(x_p(k,:,end));
                x_p_cumPeriod = zeros(N_total_period,1);
%                 for i = 1:N_total_period
%                     x_p_cumPeriod(i) = sum(x_p(k,1:i,end));
%                 end
                x_p_cumPeriod = cumsum(x_p(k,:,end));

                P_left = max(find(x_p_cumPeriod - min(data_pointer)<0));  % find the i_th period in the left of data_pointer;
                P_right = min(find(x_p_cumPeriod - max(data_pointer)>0)); % find the j_th period in the right of data_pointer;
            if (isempty(P_left)==1)+(P_right<=1)>0
                y=x_p;
                return;
            end
                
                
            if P_right+5<m2
                for i = P_right:P_right+5
                    x_p(k,i,end) = round(x_p(k,i-1,end) - seita_0(1) + 2*seita_0(1)*rand(1));
                    if (x_p(k,i,end)<40)||(x_p(k,i,end)>200)
                        x_p(k,i,1) = x_p(k,i-1,1);
                    end
                    
                    x_p(k,i,1:p_order) = x_p(k,i-1,1:p_order) + seita_0(3)*randn(1,1,p_order);
                    x_p(k,i,(p_order+1):(end-1)) = x_p(k,i-1,(p_order+1):(end-1)) + seita_0(2)*randn(1,1,K_order).*x_p(k,i-1,(p_order+1):(end-1));
                    

%                     if (x_p(k,i,2)<0.4)||(x_p(k,i,2)>1.5)
%                         x_p(k,i,2) = x_p(k,i-1,2);
%                     end       
                end
                
            else
                for i = P_right:m2
                    x_p(k,i,end) = round(x_p(k,i-1,end) - seita_0(1) + 2*seita_0(1)*rand(1));
                    if (x_p(k,i,end)<40)||(x_p(k,i,end)>200)
                        x_p(k,i,1) = x_p(k,i-1,1);
                    end
                    
                    x_p(k,i,1:p_order) = x_p(k,i-1,1:p_order) + seita_0(3)*randn(1,1,p_order);
                    x_p(k,i,(p_order+1):(end-1)) = x_p(k,i-1,(p_order+1):(end-1)) + seita_0(2)*randn(1,1,K_order).*x_p(k,i-1,(p_order+1):(end-1));
                    
              
                end
            end
        end   
        y=x_p;       
        
        
end