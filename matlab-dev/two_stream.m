kp=200;
dz=50;
z=0:dz:(kp)*dz; % z(0.5) to z(kp+0.5)
zn=z-dz/2; % z(0) to z(kp)

shortwave=true;
% shortwave=false;
albedo=0.21;
emiss=.925;
clouds=true;
cos_theta_s=0.42;

gam1=10/1000.;
gam2=0.;
tsurf=290;
g1=0.85;
tau=zeros([kp+1,1]);
t=zeros([kp+1,1]);
sigma_ray=zeros([kp+1,1]);
sigma_tot=zeros([kp+1,1]);
sigma_s_clouds=zeros([kp+1,1]);
sigma_a_clouds=zeros([kp+1,1]);
omg_s=zeros([kp+1,1]);
gamma1=zeros([kp+1,1]);
gamma2=zeros([kp+1,1]);
gamma3=zeros([kp+1,1]);
ga=zeros([kp+1,1]);
blt=zeros([kp+1,1]);

a=zeros([2*kp+1,1]);
b=zeros([2*kp+2,1]);
c=zeros([2*kp+1,1]);
r=zeros([2*kp+2,1]);
% u=zeros([2*kp+4,1]);

t(:)=tsurf-gam1.*z;

lay1=1000;
lay2=1200;
if (clouds)
    ind=find(zn>=lay1 & zn<lay2);
    sigma_s_clouds(ind)=pi./2.*100e6.*15e-6.^2.*(z(ind)-lay1)./(lay2-lay1);
%     sigma_s_clouds(ind)=pi./2.*100e6.*25e-6.^2;
    if shortwave
        sigma_a_clouds(ind)=sigma_s_clouds(ind).*0.;
    else
        sigma_a_clouds(ind)=sigma_s_clouds(ind).*.01;
    end
    sigma_s_clouds(ind)=sigma_s_clouds(ind)-sigma_a_clouds(ind);
%     t(ind(end)+1:end)=t(ind(end));
end


% set values
if shortwave
    sflux=1300.0;
%     blt(:)=(1./(exp(6.63e-34.*3e8./(10e-6.*1.381e-23.*t))-1)); 
%     blt(:)=blt(:)./blt(1).*100.;
else
    sflux=0.0;
    blt(:)=(1./(exp(6.63e-34.*3e8./(10e-6.*1.381e-23.*t))-1)); 
    blt(:)=blt(:)./blt(1).*100.;
end
sigma_ray(:)=1e-20;
sigma_tot(:)=sigma_ray(:)+sigma_s_clouds(:)+sigma_a_clouds(:);
% sigma_tot(:)=mean(sigma_tot);
omg_s(:)=(sigma_ray(:)+sigma_s_clouds(:))./(sigma_tot(:));
omg_s(:)=min(omg_s,1.0-eps*4);
ga(:)=(sigma_s_clouds(:).*g1)./(sigma_s_clouds(:)+sigma_ray(:));

% tau on half-levels: 
% zn(1)<z(1)<zn(2) => z(1)<zn(2)< z(2)=> (tau(2)-tau(1))/dz = -sigman(2)
% tau(1)=tau(2)+sigma(2)*dz
tau(kp+1)=1e-2;
for k=kp:-1:1
    tau(k)=tau(k+1)+sigma_tot(k+1)*dz;
end
% tau(1) is tau(1/2) - i.e. ground


mu1=1./sqrt(3);
% mu1=0.5;
for k=1:kp+1 % gamma(0) to gamma(kp)
   gamma1(k)=(1-omg_s(k)*(1+ga(k))*0.5)/mu1;
   gamma2(k)=omg_s(k)*(1-ga(k))/(2*mu1);
   gamma3(k)=(1-3*ga(k)*mu1*cos_theta_s)*0.5;
%    gamma1(k)=(7-omg_s(k)*(4+3*ga(k)))*0.25;
%    gamma2(k)=-(1-omg_s(k)*(4-3*ga(k)))*0.25;
%    gamma3(k)=(2-3*ga(k)*cos_theta_s)*0.25;
end


for k=2:kp+1
   Ac=(1+0.5*dz*gamma1(k)*sigma_tot(k));
   Bc=(-0.5*dz*gamma2(k)*sigma_tot(k));
   Cc=(-1+0.5*dz*gamma1(k)*sigma_tot(k));
   Dc=Bc;
   
   alp=-Bc;
   bet=-Cc;
   gam=-Bc;
   del=-Ac;
   
   Skp=0.;
   Skm=0.;
   
   if shortwave
       Skp=dz*sigma_tot(k)*gamma3(k)*omg_s(k)* ...
           sflux*exp(-0.5.*(tau(k)+tau(k-1))/cos_theta_s);
       Skm=-dz*sigma_tot(k)*(1-gamma3(k))*omg_s(k)* ...
           sflux*exp(-0.5.*(tau(k)+tau(k-1))/cos_theta_s); %+...
            %dz*sigma_tot(k)*sflux*exp(-0.5.*(tau(k)+tau(k-1))/cos_theta_s);
   end
   if ~shortwave
       Skp=Skp+dz.*sigma_tot(k)*2.*pi.*(1.-omg_s(k)).*blt(k);
       Skm=Skm-dz.*sigma_tot(k)*2.*pi.*(1.-omg_s(k)).*blt(k);
   end
   
   Ak=(Cc-Bc*gam/bet);
   Bk=(Dc-Bc*del/bet);
   Ck=(Ac-Bc*alp/bet);
   Dk=Skp-Bc/bet*Skm;

   Ek=(del-Dc*gam/Cc);
   Fk=(alp-Ac*gam/Cc);
   Gk=(bet-Bc*gam/Cc);
   Hk=Skm-gam*Skp/Cc;
   
   a(k*2-3)=Ak;
   a(2*k-2)=Ek;
   c(2*k-2)=Ck;
   c(2*k-1)=Gk;
   b(2*k-2)=Bk;
   b(2*k-1)=Fk;
   r(2*k-2)=Dk;
   r(2*k-1)=Hk;
end
b(1)=1;
c(1)=0.;
c(1)=-albedo;
r(1)=0.;
if shortwave
    r(1)=albedo*cos_theta_s*sflux*exp(-tau(1)/cos_theta_s); % tau(1) = tau(0.5)
end
if ~shortwave
    r(1)=r(1)+emiss*pi*blt(1);
end
% c(end)=1;
a(end)=-1;
b(end)=1;
r(end)=0;

% r(end)=cos_theta_s*sflux*exp(-tau(kp+1)/cos_theta_s); % tau(k+0.5)

M=diag(b)+diag(a,-1)+diag(c,1);
u=M\r;

flux_u=u(1:2:end);
flux_d=u(2:2:end);
% flux_d(1:end)=flux_d(1:end)+cos_theta_s.*sflux.*exp(-[tau]./cos_theta_s);


