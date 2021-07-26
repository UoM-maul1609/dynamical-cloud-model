function rad_tests


lambda1=500e-9;
r=logspace(log10(0.03.*lambda1./1000),log10(0.03.*lambda1),100);

clear i;
Qa=tyndall_absorber(lambda1,r,1.34-1e-9i); % needs an imaginary part of refractive index
Qs=tyndall_scatterer(lambda1,r,1.34-1e-9i); 

plot(r.*1e6,Qa);
hold on;
plot(r.*1e6,Qs,'r');
hold off;



function Q=tyndall_absorber(lambda1,r,m_lambda)

% equation 9.44, Jacobson
Q=-4.*2.*pi.*r./lambda1.*imag(((m_lambda^2-1)./(m_lambda^2+2))^2);
% Q=2.*pi.*r./lambda1.*(24*real(m_lambda)*imag(-m_lambda)/ ...
%     ((real(m_lambda)^2+imag(-m_lambda)^2)^2+4*(real(m_lambda)^2-imag(-m_lambda)^2+1)));

% Q=2.*pi.*r./lambda1.*(24*real(m_lambda)*imag(-m_lambda)/ ...
%     ((real(m_lambda)^2+2)^2));

function Q=tyndall_scatterer(lambda1,r,m_lambda)

% equation 9.46, Jacobson
Q=8/3*(2*pi.*r./lambda1).^4.*abs((m_lambda.^2-1)./(m_lambda.^2+2)).^2;
