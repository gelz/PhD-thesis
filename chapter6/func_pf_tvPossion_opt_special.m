function [lambda,err] = func_pf_tvPossion_opt_special(order_flow,Np,hyper_para)
% this function output a filtered order flow intensities given an input of
% order flow. 

    if nargin < 3
        hyper_para(1) = 1/100;
        if nargin < 2
            Np = 100;
        end
    elseif nargin < 2
        Np = 100;
    end
    
    block_size = 1;
    %hyper_para(1) controls the probablity of a jump happens for each
    %observation.
    %hyper_para(2) controls the spread of jumps...(variance)
    
    Maximum_cp = 2/hyper_para(1);  % the maximum number of change points allowed.
    hyper_para(2) = 1000;
    hyper_para(3) = 0.1;  % controls the observation noise of the true intensity;
    hyper_para(4) = 0.5;  % control the exponential decrease of order flow; 
    % observed intensity = true_intensity +
    % norm(0, hyper_para(3)*abs(true_intersitvie)); 
    
%     change_points = nan(Maximum_cp,1);
    
    [m,n] = size(order_flow);
    
    lambda = nan(m,1);
    err = lambda;
    t_train = min([500,m]);
    
    lambda(1:t_train,1) = mean(order_flow(1:t_train));
    err(1:t_train,1) = 0;
    x_p = initailize_xp(Np,lambda(1),m,t_train,Maximum_cp); %m: the length of x_p; 

    w_old = 1/Np*ones(Np,1);

    [x_p,lambda_mean,err_est] = update_xp(x_p, order_flow,t_train,w_old,hyper_para,block_size);

    lambda = lambda_mean;
    err = err_est;
    if lambda<10
        lambda = 10;
    end
end

function x_p = initailize_xp(Np,lambda,m,t_train,Maximum_cp)
    x_p = nan(Np,m);
    x_p(:,1) = normrnd(lambda,20,Np,1);
%     c_p = nan(Np,Maximum_cp);       % c_p marks the responding change points for every particle.
%     c_p(:,1) = 1;
end

function [x_p_update,lambda_mean,err_est] = update_xp(x_p,order_flow,t_train,w_old,hyper_para,block_size)
    [m,n] = size(order_flow);
    [N,p] = size(x_p);
    x_p_update = x_p;
    w_new = w_old;
    lambda_mean = nan(1,m);
    err_est = lambda_mean;
%     posterior_prob = zeros(N,m);
    log_posterior_prob = zeros(N,m);
    lambda_mean(1) = x_p(1,1);   % set it as the first value in x_p;
    err_est(1) = 0;
    current_estimate = lambda_mean(1);
    w_new = w_old;
    for i = 2:m
%      i,
%        w_new = w_old;
%         sigma_temp = fun_sample_lamda_inc(x_p_update(:,max(i-50,1):i-1),order_flow(max(i-1,1):i),hyper_para,N);
%         sigma_temp = fun_sample_lamda_inc(x_p_update(:,max(i-50,1):i-1),order_flow(max(i-5,1):i),hyper_para,N);
        sigma_temp = fun_sample_lamda_inc(x_p_update(:,max(i-50,1):i-1),order_flow(i),hyper_para,N);
       
%        for j = 1 : N 
% %            sigma_temp = fun_sample_lamda_inc(hyper_para);
%            x_p_update(j,i) = max(x_p_update(j,i-1) + sigma_temp(j),0.0001);
%            
% %            posterior_prob(j,i) = x_p_update(j,i)^order_flow(i)*exp(-x_p_update(j,i))/factorial(order_flow(i));
% %            log_posterior_prob(j,i) = -x_p_update(j,i)+order_flow(i)*log(x_p_update(j,i))-log(factorial(order_flow(i)));
%             if order_flow(i) ==0
%                 log_posterior_prob(j,i) = -x_p_update(j,i);
%             else
%                 log_posterior_prob(j,i) = -x_p_update(j,i)+order_flow(i)*log(x_p_update(j,i))-order_flow(i)*log(order_flow(i))+order_flow(i);
%             end
%        end
       
        x_p_update(:,i) = x_p_update(:,i-1) +sigma_temp;
%         x_p_update(:,i) =  x_p_update(:,i)*exp(-hyper_para(4));
        x_p_update((x_p_update(:,i)<=0),i)= 10;
        
        sigma_ob = max(hyper_para(3).*x_p_update(:,i),10);
        sigma_ob2 = max(sqrt(x_p_update(:,i)),10);
        
        new_var = 1./(1./sigma_ob.^2+1./sigma_ob2.^2);
        new_mean = (x_p_update(:,i)./sigma_ob.^2+order_flow(i)./sigma_ob2.^2).*new_var;
        
        log_posterior_prob(:,i) = 1/2.*log(new_var) - (new_mean - x_p_update(:,i)).^2./(2*sigma_ob.^2)-(new_mean-order_flow(i)).^2./(2*sigma_ob2.^2);
%         intensity_ob = (sqrt((x_p_update(:,i)-sigma_ob.^2).^2+4.*order_flow(i).*sigma_ob.^2)+(x_p_update(:,i)-sigma_ob.^2))./2;
%         if order_flow(i) ==0
%             log_posterior_prob(:,i) = -x_p_update(:,i);
%         else
%             log_posterior_prob(:,i) = -x_p_update(:,i)+order_flow(i).*log(x_p_update(:,i))-order_flow(i).*log(order_flow(i))+order_flow(i);
%         end
%         log_posterior_prob(:,i) = -intensity_ob + order_flow(i).*log(intensity_ob) - (intensity_ob - x_p_update(:,i)).^2./2./sigma_ob.^2;


%        w_new = posterior_prob(:,i).*w_old;
%        w_old = w_new; 
%        i,
            log_posterior_sum = log_posterior_prob(:,i);
            
%             log_w_new = log_posterior_sum + log(w_old);
            log_w_new = log_posterior_sum + log(w_old);
            
            log_w_new = log_w_new - max(log_w_new);
            
            w_new = exp(log_w_new);
%             max(log_posterior_sum),
%             w_new = exp(log_posterior_sum - max(log_posterior_sum));
% 
%             w_new = w_new.*w_old;
            
%             if sum(w_new ==0)
%                 w_new,
%             end
            
            w_new = w_new./sum(w_new);
        resample_decision = 0;
            if (i>500)&&(mod(i,block_size)==0)
            
                ess = 1./(sum(w_new.^2));
                if (ess > N/2)||(i<100)
                    resample_decision = 0;
                else resample_decision = 1;
                end
            end
    
        if resample_decision == 0
            w_old = w_new;
        else
            w_old = w_new;
            flags = resample([1:N]',w_old);
            for k = 1:N
                x_p_update(k,i) = x_p_update(flags(k),i);
                log_posterior_prob(k,i) = log_posterior_prob(flags(k),i);
%                 x_p(k,:,:) = x_p(flags(k),:,:);
%                 K1{k,1} = K1{flags(k),1};
%                 speech_est(k,:)=speech_est(flags(k),:);
%                 drive_est(k,:)=drive_est(flags(k),:);
%                 Pointer_num(k) = Pointer_num(flags(k)) ;
%                 Pointer_previous(k) = Pointer_previous(flags(k));
%                 Pointer_current(k) = Pointer_current(flags(k)) ;
            end
            w_old = 1./N.*ones(N,1);
        end
        current_estimate =  w_old'*x_p_update(:,i);
        err_est(i) = exp(-hyper_para(4))*realsqrt(w_old'*(x_p_update(:,i)-current_estimate).^2);
%         current_estimate = mean(mean(x_p_update(w_old ==max(w_old),i)));
        lambda_mean(i) = current_estimate*exp(-hyper_para(4));
        
       end
       
       
       


end

% function sigma_temp = fun_sample_lamda_inc(hyper_para,N)
% sigma_temp = fun_sample_lamda_inc(i,x_p_update(:,max(i-50,1):i),order_flow(max(i-50,1):i),hyper_para,N);

function        sigma_temp = fun_sample_lamda_inc(x_p_update,order_flow,hyper_para,N)
    k = rand(N,1);
    index = find(k>1-hyper_para(1));
    index2 = find(k<=1-hyper_para(1));
    sigma_temp  = zeros(N,1);
    
    sigma_temp(index) = mean(order_flow)-x_p_update(index,end)+ max(sqrt(abs(order_flow(end))),10)*randn(length(index),1);
%     sigma_temp(index2) = x_p_update(index2,end).*(exp(-hyper_para(4))-1);
    
%     sigma_temp = hyper_para(2)*(rand(N,1)-0.5);
%     sigma_temp(k<1-hyper_para(1)) =0;
    
    
    
%     if k<1-hyper_para(1)
%         sigma_temp = 0;
%     else
%     end

end




