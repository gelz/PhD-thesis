  function [p_inc] = alg_pdf_fftv7(x,y,paras)
% add speed up for large x and y's.
% use power and polar form to see if it can be speed up with better
% accuracy;

% if nargin <1
%     x = 1;
%     y = 1;
%     paras = [1.0 1.2 1.0 1.2];  % paras = [p01 p0-1 p10 p-10]
% elseif nargin < 3
%     paras = [1.0 1.2 1.0 1.2];
% end

p_01 = paras(1)./sum(paras);
p_02 = paras(2)./sum(paras);  % p_02 denote for p_0-1
p_10 = paras(3)./sum(paras);
p_20 = paras(4)./sum(paras);  % p_20 denote for p_-10

% Z = quadgk(@cce_CONT_RW_term2,0,pi);  

    function  z = cce_Characteristic(f)
        y1 = character_function(p_01,p_02,-f);
        y2 = character_function(p_10,p_20,f);
%         z = y1.^x.*y2.^y./1i./f;  % this is right: with 2*pi scalings of T_shift;
        z1 = y1.^x.*y2.^y;
%         z1 = r1.^x.*r2.^y.*sin(x.*theta1+y.*theta2);
        z = imag(z1)./f;
    end
% t = -1:0.00001:1;
% z = cce_Characteristic(t);
% figure();plot(t,real(z));

% p = 0;
if max(x,y) < 20
    bound_max = 100;
else 
    bound_max = 50;
end
%     p_inc = 0.5 - 1./2/pi.*(real(quadgk(@cce_Characteristic,-pi/2,-0.00001))+real(quadgk(@cce_Characteristic,0.00001,pi/2)));
    p_inc = 0.5 - 1./pi.*(quadgk(@cce_Characteristic,0,bound_max));
%     p_inc = 0.5 - 1./pi.*(integral(@cce_Characteristic,-bound_max,0));



end




function y = character_function(pxy,pxy2,u)

a = (pxy2+pxy)./(pxy+pxy2-1i*u)*(pxy)./(pxy+pxy2);
c = (pxy2+pxy)./(pxy+pxy2-1i*u)*(pxy2)./(pxy+pxy2);

y = (1-sqrt(1-4.*a.*c))./2./a;

% r = abs(y);
% theta = angle(y);

end


