% See Toon et al. pp16292
N=100; % vertical dimension of atmosphere
dz=50;
M=2*N; % size of matrix to solve

cloudb=2000;
cloudt=2200;
cloudn=100e6;
cloudd=15e-6;

cos_theta_s=0.4;
albedo=0.2;
g1=0.85;
sflux=1365;

model='eddington';
% model='quadrature';
delta=true;


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

% tri-diagonal terms
a=zeros([M-1,1]);
b=zeros([M,1]);
c=zeros([M-1,1]);
r=zeros([M,1]);


% rayleigh scattering - small
sigma_ray(:)=1e-10+z*1e-10;

% calculate the cloud extinction
for i=1:N-1
    if((z(i+1)>=cloudb) && (z(i)<cloudt))
        sigma_s_cloud(i)=pi/2*cloudn*cloudd^2;
    end
end

% absorption
sigma_a_cloud(:)=sigma_s_cloud(:).*0.0001;
sigma_s_cloud(:)=sigma_s_cloud(:)-sigma_a_cloud(:);

% calculate total extinction
sigma_tot(:)=sigma_ray(:)+sigma_s_cloud(:)+sigma_a_cloud(:);
% calculate the single scattering albedo
omg_s(:)=(sigma_ray(:)+sigma_s_cloud(:))./sigma_tot(:);
omg_s(:)=min(omg_s(:),1-10.*eps);
% asymmetry factor
ga(:)=(sigma_s_cloud(:).*g1)./(sigma_s_cloud(:)+sigma_ray(:));

% calculate optical depth in each layer
% plot(flipud(cumsum(flipud(tau))),z)
tau(:)=sigma_tot(:).*dz;
% cumululative optical depth above layer n
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

% calculate lambdas - equation 21
for i=1:N
    lambda1(i)=(gamma1(i)^2-gamma2(i)^2)^0.5;
end
% calculate Gammas - equation 22
for i=1:N
%     Gamma1(i)=gamma2(i)./(lambda1(i)+gamma1(i));
    Gamma1(i)=(gamma1(i)-lambda1(i))./(gamma2(i));
end

% set the incoming diffusive flux to zero
Fm00=0;
% equation 24
Cm10=omg_s(1).*sflux.*exp(-(tauc(1))./cos_theta_s).* ...
    ((gamma1(1)+1./cos_theta_s).*gamma4(1)+gamma3(1)*gamma2(1))/ ...
    (lambda1(1)^2-1./cos_theta_s^2);

% first fill the "solution" using equation 41
for i=1:M
    if (mod(i,2)==0) % even
        if(i<=M-2) 
            n=i/2; % equation 40?

            a(i)=(1-Gamma1(n+1)*exp(-lambda1(n+1)*tau(n+1))).* ...
                (1-Gamma1(n)*exp(-lambda1(n)*tau(n))) - ...
                (Gamma1(n)+exp(-lambda1(n)*tau(n))).* ...
                (Gamma1(n+1)-exp(-lambda1(n+1)*tau(n+1)));
            c(i)=(1+Gamma1(n+1)*exp(-lambda1(n+1)*tau(n+1))).* ...
                (Gamma1(n+1)-exp(-lambda1(n+1)*tau(n+1))) - ...
                (1-exp(-lambda1(n+1)*tau(n+1))).* ...
                (1+exp(-lambda1(n+1)*tau(n+1)));
            b(i)=(1-Gamma1(n)*exp(-lambda1(n)*tau(n))).* ...
                (1-Gamma1(n+1)*exp(-lambda1(n+1)*tau(n+1))) - ...
                (Gamma1(n)-exp(-lambda1(n)*tau(n))).* ...
                (Gamma1(n+1)-exp(-lambda1(n+1)*tau(n+1)));
            
            Cpnp10=omg_s(n+1)*sflux* ...
                exp(-(tauc(n+1)+0)/cos_theta_s)*...
                ((gamma1(n+1)-1/cos_theta_s)*gamma3(n+1)+...
                gamma4(n+1)*gamma2(n+1)) ...
                /(lambda1(n+1)^2-1./cos_theta_s^2);
            Cpnt=omg_s(n)*sflux* ...
                exp(-(tauc(n)+tau(n))/cos_theta_s)*...
                ((gamma1(n)-1/cos_theta_s)*gamma3(n)+...
                gamma4(n)*gamma2(n)) ...
                /(lambda1(n)^2-1./cos_theta_s^2);
            Cmnp10=omg_s(n+1)*sflux* ...
                exp(-(tauc(n+1)+0)/cos_theta_s)*...
                ((gamma1(n+1)+1/cos_theta_s)*gamma4(n+1)+...
                gamma3(n+1)*gamma2(n+1)) ...
                /(lambda1(n+1)^2-1./cos_theta_s^2);
            Cmnt=omg_s(n)*sflux* ...
                exp(-(tauc(n)+tau(n))/cos_theta_s)*...
                ((gamma1(n)+1/cos_theta_s)*gamma4(n)+...
                gamma3(n)*gamma2(n)) ...
                /(lambda1(n)^2-1./cos_theta_s^2);
            
            % e2n+1*(C+n+1(0) - Cn+(tau))+e4n+1*(C-n+1(0)-C-n(tau))
            r(i)=(1-Gamma1(n)*exp(-lambda1(n)*tau(n))).* ...
                (Cpnp10-Cpnt) + ...
                (Gamma1(n+1)-exp(-lambda1(n+1)*tau(n+1))).* ...
                (Cmnp10-Cmnt);
        end
    else % odd
        % note, I think it should be up to M-2 inclusive, i.e. this
        if(i<=M) 
        
            n=(i+1)/2; % equation 40?

            if i==1
                a(i)=0.;
                b(i)=1+gamma1(n)*exp(-lambda1(n)*tau(n));
                c(i)=-(1-gamma1(n)*exp(-lambda1(n)*tau(n)));
                r(i)=Fm00-Cm10;
            else
                a(i)=(1-Gamma1(n)*exp(-lambda1(n)*tau(n))).* ...
                    (Gamma1(n)+exp(-lambda1(n)*tau(n))) - ...
                    (Gamma1(n)-exp(-lambda1(n)*tau(n))).* ...
                    (1+Gamma1(n)*exp(-lambda1(n)*tau(n)));
                if(n<N) 
                    c(i)=(Gamma1(n)+exp(-lambda1(n)*tau(n))).* ...
                        (Gamma1(n+1)-exp(-lambda1(n+1)*tau(n+1))) - ...
                        (1+exp(-lambda1(n)*tau(n))).* ...
                        (1-exp(-lambda1(n+1)*tau(n+1)));
                    % e1n*e1,n+1 - e3n*e3,n+1
                    b(i)=(1+Gamma1(n)*exp(-lambda1(n)*tau(n))).* ...
                        (1+Gamma1(n+1)*exp(-lambda1(n+1)*tau(n+1))) - ...
                        (Gamma1(n)+exp(-lambda1(n)*tau(n))).* ...
                        (Gamma1(n+1)+exp(-lambda1(n+1)*tau(n+1)));
                    Cpnp10=omg_s(n+1)*sflux* ...
                        exp(-(tauc(n+1)+0)/cos_theta_s)*...
                        ((gamma1(n+1)-1/cos_theta_s)*gamma3(n+1)+...
                        gamma4(n+1)*gamma2(n+1)) ...
                        /(lambda1(n+1)^2-1./cos_theta_s^2);
                    Cpnt=omg_s(n)*sflux* ...
                        exp(-(tauc(n)+tau(n))/cos_theta_s)*...
                        ((gamma1(n)-1/cos_theta_s)*gamma3(n)+...
                        gamma4(n)*gamma2(n)) ...
                        /(lambda1(n)^2-1./cos_theta_s^2);
                    Cmnp10=omg_s(n+1)*sflux* ...
                        exp(-(tauc(n+1)+0)/cos_theta_s)*...
                        ((gamma1(n+1)+1/cos_theta_s)*gamma4(n+1)+...
                        gamma3(n+1)*gamma2(n+1)) ...
                        /(lambda1(n+1)^2-1./cos_theta_s^2);
                    Cmnt=omg_s(n)*sflux* ...
                        exp(-(tauc(n)+tau(n))/cos_theta_s)*...
                        ((gamma1(n)+1/cos_theta_s)*gamma4(n)+...
                        gamma3(n)*gamma2(n)) ...
                        /(lambda1(n)^2-1./cos_theta_s^2);

                    % e3n*(C+n+1(0) - Cn+(tau))+e1n*(C-n(tau)-C-n+1(0))
                    r(i)=(Gamma1(n)+exp(-lambda1(n)*tau(n))).* ...
                        (Cpnp10-Cpnt) + ...
                        (1+Gamma1(n)*exp(-lambda1(n)*tau(n))).* ...
                        (Cmnt-Cmnp10);
                end
            end
        end
    end
    if(i==2*N)
        % equation 43
        n=i/2;
        c(end)=0.;
        a(end)=(1+Gamma1(n)*exp(-lambda1(n)*tau(n)))-...
            albedo*(Gamma1(n)+exp(-lambda1(n)*tau(n)));
        % e2N - Rsfc*e4N
        b(end)=1-Gamma1(n)*exp(-lambda1(n)*tau(n)) - ...
            albedo*(Gamma1(n)-exp(-lambda1(n)*tau(n)));

        CpNt=omg_s(n)*sflux* ...
            exp(-(tauc(n)+tau(n))/cos_theta_s)*...
            ((gamma1(n)-1/cos_theta_s)*gamma3(n)+...
            gamma4(n)*gamma2(n)) ...
            /(lambda1(n)^2-1./cos_theta_s^2);
        CmNt=omg_s(n)*sflux* ...
            exp(-(tauc(n)+tau(n))/cos_theta_s)*...
            ((gamma1(n)+1/cos_theta_s)*gamma4(n)+...
            gamma3(n)*gamma2(n)) ...
            /(lambda1(n)^2-1./cos_theta_s^2);
        
        r(end)=-CpNt+albedo*CmNt;
    end
end

M1=diag(b)+diag(a,-1)+diag(c,1);
u=M1\r;





