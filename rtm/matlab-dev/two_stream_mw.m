% See Toon et al. pp16292
N=1000; % vertical dimension of atmosphere
dz=10;
M=2*N+2; % size of matrix to solve

shortwave=true;
shortwave=false;
model='eddington';
% model='quadrature';
delta=true;
% delta=false;

% cloudb=2000;
% cloudt=7001;
cloudb=2000;
cloudt=2101;
cloudn=100e6;
cloudd=20e-6;
blt=100;

cos_theta_s=0.40;
albedo=0.2;
emiss=0.925;
g1=0.85;
sflux=1365;
tsurf=290;
gam1=10./1000;

% model variables
z=(0:dz:(N-1)*dz)';
z=flipud(z);
tau=zeros([N,1]);
sigma_ray=zeros([N,1]);
sigma_tot=zeros([N,1]);
sigma_s_cloud=zeros([N,1]);
sigma_a_cloud=zeros([N,1]);
omg_s=zeros([N,1]);
ga=zeros([N,1]);
gamma1=zeros([N,1]);
gamma2=zeros([N,1]);
gamma3=zeros([N,1]);
gamma4=zeros([N,1]);
lambda1=zeros([N,1]);
Gamma1=zeros([N,1]);
t=zeros([N,1]);


% tri-diagonal terms
a=zeros([M-1,1]);
b=zeros([M,1]);
c=zeros([M-1,1]);
r=zeros([M,1]);

% temperature
t(:)=tsurf-gam1.*z;

% rayleigh scattering - small
sigma_ray(:)=1e-20;

% calculate the cloud extinction
for i=1:N-1
    if((z(i+1)>=cloudb) && (z(i)<cloudt))
        sigma_s_cloud(i)=pi/2*cloudn*cloudd^2.*(z(i)-cloudb)./(cloudt-cloudb);
        sigma_s_cloud(i)=pi/6*1000*cloudn*cloudd^3.*(z(i)-cloudb)./(cloudt-cloudb);
        sigma_s_cloud(i)=cloudn.*((sigma_s_cloud(i)/cloudn).*6./(pi.*1000)).^(2/3);
%         sigma_s_cloud(i)=0.0040/2.;
    end
end

% absorption
if shortwave
    sigma_a_cloud(:)=sigma_s_cloud(:).*0.;
else
    sigma_a_cloud(:)=sigma_s_cloud(:).*0.001;    
end
sigma_s_cloud(:)=sigma_s_cloud(:)-sigma_a_cloud(:);

% calculate total extinction
sigma_tot(:)=sigma_ray(:)+sigma_s_cloud(:)+sigma_a_cloud(:);
% calculate the single scattering albedo
omg_s(:)=(sigma_ray(:)+sigma_s_cloud(:))./sigma_tot(:);
omg_s(:)=min(omg_s(:),1-4.*eps);
% asymmetry factor
ga(:)=(sigma_s_cloud(:).*g1)./(sigma_s_cloud(:)+sigma_ray(:));

% calculate optical depth in each layer
% plot(flipud(cumsum(flipud(tau))),z)
tau(:)=sigma_tot(:).*dz;
tauc=cumsum(tau)-tau;

if delta
    % delta=scaling
    for i=1:N
        tau(i)=(1-omg_s(i)*ga(i).^2).*tau(i);
        omg_s(i)=(1-ga(i)^2)*omg_s(i)/(1-omg_s(i)*ga(i)^2);
        ga(i)=ga(i)./(1+ga(i));
    end
    % cumululative optical depth above layer n
    tauc=cumsum(tau)-tau;
end


switch model
    case 'quadrature'
        mu1=1./sqrt(3);
        for i=1:N % parameters in equations
           gamma1(i)=(1-omg_s(i)*(1+ga(i))*0.5)/mu1;
           gamma2(i)=omg_s(i)*(1-ga(i))/(2*mu1);
           gamma3(i)=(1-3*ga(i)*mu1*cos_theta_s)*0.5;
           gamma4(i)=1.-gamma3(i);
        end        
    case 'eddington'
        mu1=0.5;
        for i=1:N % parameters in equations
           gamma1(i)=(7-omg_s(i)*(4+3*ga(i)))*0.25;
           gamma2(i)=-(1-omg_s(i)*(4-3*ga(i)))*0.25;
           gamma3(i)=(2-3*ga(i)*mu1)*0.25;
           gamma4(i)=1.-gamma3(i);
        end
        
    otherwise
        disp('unknown');
end

if ~shortwave
    gamma3(:)=0.;gamma4(:)=0;
end


% calculate reflectances
for i=1:N
    k=(gamma1(i)^2-gamma2(i)^2)^0.5; % eq 18
    alpha1=gamma1(i)*gamma4(i)+gamma2(i)*gamma4(i); % eq 16
    alpha2=gamma1(i)*gamma3(i)+gamma2(i)*gamma4(i); % eq 17
    % eq 14
    b((i)*2)=-omg_s(i)/((1-k^2*cos_theta_s^2)* ...
        (k+gamma1(i))*exp(k*tau(i))+(k-gamma1(i))*exp(-k*tau(i))) * ...
        ((1-k*cos_theta_s)*(alpha2+k*gamma3(i))*exp(k*tau(i)) - ...
        (1+k*cos_theta_s)*(alpha2-k*gamma3(i))*exp(-k*tau(i)) - ...
        2*k*(gamma3(i)-alpha2*cos_theta_s)*exp(-tau(i)/cos_theta_s));
    b((i)*2+1)=b((i)*2);
    
    % eq 25
%     b((i)*2)=(gamma2(i)*(1-exp(-2*k*tau(i))))/ ...
%         (k+gamma1(i)+(k-gamma1(i))*exp(-2*k*tau(i)));
%     b((i)*2+1)=b((i)*2);
    
end
b(end)=-albedo;
% calculate transmissions
for i=1:N
    k=(gamma1(i)^2-gamma2(i)^2)^0.5; % eq 18
    alpha1=gamma1(i)*gamma4(i)+gamma2(i)*gamma4(i); % eq 16
    alpha2=gamma1(i)*gamma3(i)+gamma2(i)*gamma4(i); % eq 17
    % transmissions, eq 15
    a(i*2)=-exp(-tau(i)/cos_theta_s)* ...
        (1-omg_s(i)/((1-k^2*cos_theta_s^2)* ...
        (k+gamma1(i))*exp(k*tau(i))+(k-gamma1(i))*exp(-k*tau(i))) * ...
        ((1+k*cos_theta_s)*(alpha1+k*gamma4(i))*exp(k*tau(i)) - ...
        (1-k*cos_theta_s)*(alpha1-k*gamma4(i))*exp(-k*tau(i)) - ...
        2*k*(gamma4(i)+alpha1*cos_theta_s)*exp(tau(i)/cos_theta_s)));
    c(i*2)=a(i*2);

    % eq 26
%     a(i*2)=2*k*exp(-k*tau(i))/(k+gamma1(i)+(k-gamma1(i))*exp(-2*k*tau(i)));
%     c(i*2)=a(i*2);
% 
%     
    a(i*2-1)=1;
    c(i*2-1)=1;
end
a(end)=1;
c(end)=1;

% calculate source terms = i=1 is top
r(1)=0.;
if shortwave
    r(1)=sflux*cos_theta_s;
    r(end)=albedo*cos_theta_s*sflux*exp(-tauc(end)/cos_theta_s);
else
    b(1)=-1;
    r(1)=0 ; 
    b(end)=0.;
    r(end)=emiss*pi*blt;
end
for i=1:N
    % upward
%     r(i*2)=-sflux*exp(-(tauc(i))/cos_theta_s).*gamma3(i)*omg_s(i)*tau(i);
%     % downward
%     r(i*2+1)=sflux*exp(-(tauc(i))/cos_theta_s).*gamma4(i)*omg_s(i)*tau(i);

    if ~shortwave
        r(2*i)=2*pi*(1-omg_s(i)).*blt.*tau(i).*(1-z(end-1)./(z(i)+z(end-1)));
        r(2*i+1)=2*pi*(1-omg_s(i)).*blt*tau(i).*(1-z(end-1)./(z(i)+z(end-1)));
    end
end






% solve
M1=diag(b)+diag(a,-1)+diag(c,1);
u=M1\r;


flux_u=u(1:2:end);
flux_d=u(2:2:end);
% flux_d(1:end)=flux_d(1:end)+cos_theta_s*sflux*exp(-[0; tauc(:); ]/cos_theta_s);

if shortwave
    % see Liou, 3.5.3 (and 3.5.1a), page 108
    plot(-diff(flipud(flux_d(1:end))-flipud(flux_u(1:end)))./dz.*86400./1000,...
        flipud(z),'linewidth',2)
else
    % the IR cooling rate is the opposite
    % see Liou, 4.7.1, and 4.7.2, page 160
    plot(-diff(flipud(flux_u(1:end))-flipud(flux_d(1:end)))./dz.*86400./1000,...
        flipud(z),'linewidth',2)
end

