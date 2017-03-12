function [x_p_est,R] = joint_est_para_sin_v2(p_order,K_order,data_input,fs,T0)
    
%     x  =  data_input(2*T0+1:T0);
    x  =  data_input((1+T0):3*T0);
    unit_freq = 2*pi/T0;
    unit_phase = unit_freq *(1:length(x));
    
    for k = 1 : (K_order/2)
        c(:,k) = cos((k-1)*unit_phase);
    end
    for k = ((K_order/2)+1) : K_order
        c(:,k) = sin((k-K_order/2)*unit_phase);
    end
%     plot(c);
%     R1_tmp = corrmtx(x,p_order);
%     R1 = R1_tmp(2:end,2:end);
    
    for i = 1:p_order
        for j = 1:p_order
            R1(i,j) = cir_cov(x,x,j,i); 
        end
    end
    
    for i = 1:p_order
        for j = 1:K_order
            R2(i,j) = cir_cov(c(:,j),x,0,i);
        end
    end
    
%     for i = 1:K_order
%         R3_diag(i) = sum(c(:,i).^2);
%     end
%     
%     R3 = diag(R3_diag);
    for i = 1:K_order
        for j = 1:K_order
        R3(i,j) = cir_cov(c(:,j),c(:,i),0,0);
        end
    end
    
    for i = 1:p_order
        p(i) = -cir_cov(x,x,0,i);
    end
    for i = p_order+1:(p_order+K_order)
        p(i) = cir_cov(c(:,i-p_order),x,0,0);        
    end
    
    R = [R1,-R2; -R2', R3];
%     R,
%     x_p_est = inv(R)*p';
    x_p_est = R\p';
    
%     (x_p_est'*R*x_p_est - 2*x_p_est'*p')./length(data_input),
    
end



function x = cir_cov(x1,x2,lag1,lag2)
    x = sum(circshift(x1,lag1).*circshift(x2,lag2));
end




    