function [x,sum1]=inline_planck_curve(lambda,temp)

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
x=h_planck.*c_light./(k_boltz.*temp.*lambda);
sum1=0.;
for j=1:min(512 , 2+round(20./x)) % criteria on webpage
    sum1=sum1+(x.^3./real(j) + ...
            3.*x.^2./real(j).^2 + ...
            6.*x./real(j).^3 + ...
            6./real(j).^4) .* ...
                exp(-real(j).*x);
end



