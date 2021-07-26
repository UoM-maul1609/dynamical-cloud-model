function [blt,x_upper,x_lower,sum_upper,sum_lower]=inline_planck(lambda_low,lambda_high,temp)

c_light=3e8;
k_boltz=1.381e-23;
h_planck=6.63e-34;
fac_planck=2.*k_boltz.^4./(h_planck.^3.*c_light.^2).*pi;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% calculate the planck function in these bins				 				 !
% see http://www.spectralcalc.com/blackbody/inband_radiance.html             !
% note, the integral result is the same no matter form of BB function        !
% but the value of x is what is in the exponential                           !
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fac=fac_planck .* temp.^4;
x_upper=h_planck.*c_light/(k_boltz.*temp.*lambda_high);
x_lower=h_planck.*c_light/(k_boltz.*temp.*lambda_low);
sum_upper=0.;
for j=1:min(512 , 2+round(20./x_upper)) % criteria on webpage
    sum_upper=sum_upper+(x_upper.^3./real(j) + ...
            3.*x_upper.^2./real(j).^2 + ...
            6.*x_upper/real(j).^3 + ...
            6./real(j).^4) * ...
                exp(-real(j).*x_upper);
end

sum_lower=0.;
for j=1:min(512 , 2+round(20./x_lower)) % criteria on webpage
    sum_lower=sum_lower+(x_lower.^3./real(j) + ...
            3.*x_lower.^2./real(j).^2 + ...
            6.*x_lower/real(j).^3 + ...
            6./real(j).^4) * ...
                exp(-real(j).*x_lower);
end
blt = fac.*(sum_upper-sum_lower);

