function [T_out2,std_T_out2] = mean_out_TimeVary_harmonics_RB(w_new,x_p,i);   


        [m1,m2,m3] = size(x_p);
        ss =0;
        ss2 = zeros(m1,1);
        std_T_out2 = 0;
        
        Time_index = i;
        for k=1:m1
                
                Period_num = length(find(x_p(k,:,end)~=0));
        
                
                ss = w_new(k).*x_p(k,Period_num,end)+ss;
                ss2(k) = x_p(k,Period_num,end);

        end
        
        T_out2 = ss;
        std_T_out2 = sqrt(dot((ss2-ss).^2,w_new));
end