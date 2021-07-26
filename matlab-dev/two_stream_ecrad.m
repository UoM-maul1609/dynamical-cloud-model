% close all

N=500; % vertical dimension of atmosphere
dz=10;
M=2*N+2; % size of matrix to solve

shortwave=true;
% shortwave=false;
model='eddington';
% model='quadrature';
delta=true;
% delta=false;

LwDiffusivity= 1.66;
SwDiffusivity= 2.00;

% cloudb=2000;
% cloudt=7001;
cloudb=2000;
cloudt=2201;
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
zl=z;
z=z+0.5*dz;
zu=z+0.5*dz;
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
blt_top=zeros([N,1]);
blt_bot=zeros([N,1]);
source_up=zeros([N,1]);
source_dn=zeros([N,1]);
R=zeros([N,1]);
T=zeros([N,1]);

trans_dir_dir=zeros([N,1]);
ref_diff=zeros([N,1]);
trans_diff=zeros([N,1]);
ref_dir=zeros([N,1]);
trans_dir_diff=zeros([N,1]);

% tri-diagonal terms
a=zeros([M-1,1]);
b=zeros([M,1]);
c=zeros([M-1,1]);
r=zeros([M,1]);

% temperature
tl=tsurf-gam1.*zl;
t=tl-dz*gam1;
tu=t-dz*gam1;

% planck
blt_top(:)=pi.*blt.*tu.^4./(t(end).^4);
blt_bot(:)=pi.*blt.*tl.^4./(t(end).^4);

% rayleigh scattering - small
sigma_ray(:)=1e-20;

% calculate the cloud extinction
for i=1:N
    if((z(i)>=cloudb) && (z(i)<cloudt))
        sigma_s_cloud(i)=pi/2*cloudn*cloudd^2.*(z(i)-cloudb)./(cloudt-cloudb);
        sigma_s_cloud(i)=pi/6*1000*cloudn*cloudd^3.*(z(i)-cloudb)./(cloudt-cloudb);
        sigma_s_cloud(i)=cloudn.*((sigma_s_cloud(i)/cloudn).*6./(pi.*1000)).^(2/3);
    end
end

% absorption
if shortwave
    sigma_a_cloud(:)=sigma_s_cloud(:).*0.0136;
else
    sigma_a_cloud(:)=sigma_s_cloud(:).*0.481;    
end
sigma_s_cloud(:)=sigma_s_cloud(:)-sigma_a_cloud(:);

% calculate total extinction
sigma_tot(:)=sigma_ray(:)+sigma_s_cloud(:)+sigma_a_cloud(:);
% calculate the single scattering albedo
omg_s(:)=(sigma_ray(:)+sigma_s_cloud(:))./sigma_tot(:);
omg_s(:)=min(omg_s(:),1-10.*eps);
% asymmetry factor
ga(:)=(sigma_s_cloud(:).*g1)./(sigma_s_cloud(:)+sigma_ray(:));
% omg_s(:)=1;
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


if shortwave
    % PIFM https://ntrs.nasa.gov/api/citations/19860012551/downloads/19860012551.pdf
    for i=1:N
        factor=0.75*ga(i);
        gamma1(i)=2-omg_s(i)*(1.25+factor);
        gamma2(i)=omg_s(i)*(0.75-factor);
        gamma3(i)=0.5-cos_theta_s*factor;
        gamma4(i)=1-gamma3(i);
        
    end
else
    % calc_two_stream_gammas_lw - looks like quadrature method
    % but Lw is 1.66 instead of 1.73
    for i=1:N
        factor=(LwDiffusivity*0.5)*omg_s(i);
        gamma1(i)=LwDiffusivity-factor*(1+ga(i));
        gamma2(i)=factor*(1-ga(i));
    end    
end




if ~ shortwave
    for i=1:N
        if tau(i)>1.e-3
            k_exponent = sqrt(max((gamma1(i)-gamma2(i))*(gamma1(i)+gamma2(i)), ...
                1.e-12)); % equation 18 of meador and weaver (1980)
            exponential = exp(-k_exponent*tau(i));
            exponential2 = exponential*exponential;
            reftrans_factor = 1.0 / ...
                (k_exponent+gamma1(i)+(k_exponent - gamma1(i))*exponential2);
            % Meador and Weaver (1980) eq 25.
            R(i) = gamma2(i)*(1-exponential2)*reftrans_factor;
            % Meador and Weaver (1980) Eq. 26
            T(i)=2.*k_exponent*exponential*reftrans_factor;

            % upward and downward emission assuming Planck function varies lineary
            % with optical depth within the layer (Wiscombe, JQSRT 1976)

            % from Stackhouse and Stephens (JAS 1991) Eqs 5 & 12
            coeff = (blt_bot(i)-blt_top(i)) / (tau(i)*(gamma1(i)+gamma2(i)));
            coeff_up_top = coeff + blt_top(i);
            coeff_up_bot = coeff + blt_bot(i);
            coeff_dn_top =-coeff + blt_top(i);
            coeff_dn_bot =-coeff + blt_bot(i);
            source_up(i) = coeff_up_top - R(i)*coeff_dn_top - ...
                T(i)*coeff_up_bot;
            source_dn(i) = coeff_dn_bot - R(i)*coeff_up_bot - ...
                T(i)*coeff_dn_top;
        else
            % Eq. 18 of Meador and Weaver (1980)
            k_exponent = sqrt(max((gamma1(i)-...
                gamma2(i))*(gamma1(i)+gamma2(i)), 1e-12));
            R(i)=gamma2(i)*tau(i);
            T(i)=(1-k_exponent*tau(i))/(1+tau(i)*(gamma1(i)-k_exponent));
%             T(i)=1-R(i)-tau(i)/cos_theta_s*(1-omg_s(i)*(gamma3(i)+gamma4(i)));
            source_up(i)=(1-R(i)-T(i)) * ...
                0.5 * (blt_top(i)+blt_bot(i));
            source_dn(i)=source_up(i);

        end
    end
else
    % short wave radiation
    for i=1:N
        
        alpha1 = gamma1(i)*gamma4(i)    +gamma2(i)*gamma3(i); % eq. 16
        alpha2 = gamma1(i)*gamma3(i)    +gamma2(i)*gamma4(i); % eq. 17
        k_exponent = sqrt(max((gamma1(i)-gamma2(i))*(gamma1(i)+gamma2(i)), ...
            1.e-12)); % equation 18 of meador and weaver (1980)
        
        mu0_local = cos_theta_s;
        % if K*cos_theta_s is very close to 1 then ref_dif and
        % trans_dir_diff can be outside 0-1
        if(abs(1-k_exponent*cos_theta_s) < 1000*eps)
            mu0_local=cos_theta_s*(1-10*eps);
        end
        
        od_over_mu0 = max(tau(i)/mu0_local, 0);
        k_mu0 = k_exponent*mu0_local;
        k_gamma3 = k_exponent*gamma3(i);
        k_gamma4 = k_exponent*gamma4(i);
        
        % check for mu0 <= 0
        exponential0 = exp(-od_over_mu0);
        trans_dir_dir(i) = exponential0; % transmission of direct radiation to direct radiation
        exponential = exp(-k_exponent*tau(i));
        
        exponential2 = exponential*exponential;
        k_2_exponential = 2.*k_exponent*exponential;
        
        reftrans_factor = 1.0 / (k_exponent+gamma1(i) + (k_exponent - gamma1(i))*exponential2);
        
        % Meador and Weaver (1980) Eq. 25
        ref_diff(i) = gamma2(i)*(1-exponential2)*reftrans_factor;
        
        % Meador and Weaver (1980) Eq. 26
        trans_diff(i) = k_2_exponential * reftrans_factor;
        
        % Because we assume the incoming "direct" flux comes into a plane
        % perpendicular to the sun, we need to multiply by mu0 as diffuse
        % fluxes are always perpendicular to ground - Robin Hogan (Personal
        % Communication)
        reftrans_factor = mu0_local*omg_s(i)*reftrans_factor / (1-k_mu0*k_mu0);
        
        % Meador and Weaver (1980) Eq. 14. multiplying top and bottom by
        % exp(-k_exponent*tauc(i)) in case of very high opticla depths
        % reflectance to direct radiation
        ref_dir(i) = reftrans_factor *((1-k_mu0)*(alpha2+k_gamma3) - ...
            (1+k_mu0) * (alpha2 - k_gamma3)*exponential2 - ...
            k_2_exponential*(gamma3(i) - alpha2*mu0_local)*exponential0);
        
        % Meador and Weaver (1980) Eq. 15 multiplying top and bottom by
        % exp(-k_exponent*tauc(i)), minus the 1*exp(-tauc(i)/mu0) term
        % representing direct unscattered transmittance.
        trans_dir_diff(i) = reftrans_factor*(k_2_exponential* ...
            (gamma4(i)+alpha1*mu0_local) - exponential0*((1+k_mu0)* ...
            (alpha1+k_gamma4)-(1-k_mu0)*(alpha1-k_gamma4)*exponential2));
        % final check that ref_dir + trans_dir_diff <= 1
        ref_dir(i) = max(0,min(ref_dir(i),1));
        trans_dir_diff(i) = max(0,min(trans_dir_diff(i), 1-ref_dir(i)));
    end
    
end

% ecrad 
% - calc_two_stream_gammas_lw - done
% - calc_two_stream_gammas_sw - done
% - calc_reflectance_transmittance_lw - done
% -     calc_reflectance_transmittance_isothermal_lw
% -     calc_no_scattering_transmittance_lw
% - calc_reflectance_transmittance_sw
% -     calc_reflectance_transmittance_z_sw
% -     calc_frac_scattered_diffuse_sw


% if ~shortwave
%     % build the matrix - a
%     for i=1:N
%         a(i*2-1)=1;
%         a(i*2)=-T(i);
%     end
%     a(end)=1;
%     % build the matrix - c
%     for i=1:N
%         c(i*2-1)=1;
%         c(i*2)=-T(i);
%     end
%     c(end)=1;
%     % build the matrix - b
%     for i=1:N
%         b(i*2)=-R(i);
%         b(i*2+1)=-R(i);
%     end
%     b(end)=-albedo;
%     % build the matrix - r
%     for i=1:N
%         r(i*2)=source_up(i);
%         r(i*2+1)=source_dn(i);
%     end
%     r(end)=emiss*pi*blt.*tl(end).^4./(t(end).^4);;
%     
%     % solve++++++++++++++++++
%     M1=diag(b)+diag(a,-1)+diag(c,1);
%     u=M1\r;
% 
% 
%     flux_up=u(1:2:end);
%     flux_dn=u(2:2:end);
% 
% 
% end
% 
% if shortwave 
%     % direct radiation
%     R=ref_dir;
%     T=trans_dir_dir; %+trans_dir_diff;
%     % build the matrix - a
%     for i=1:N
%         a(i*2-1)=1;
%         a(i*2)=-T(i);
%     end
%     a(end)=1;
%     % build the matrix - c
%     for i=1:N
%         c(i*2-1)=1;
%         c(i*2)=-T(i);
%     end
%     c(end)=1;
%     % build the matrix - b
%     for i=1:N
%         b(i*2)=-R(i);
%         b(i*2+1)=-R(i);
%     end
%     b(end)=-albedo;
%     % build the matrix - r
%     for i=1:N
%         r(i*2)=source_up(i);
%         r(i*2+1)=source_dn(i);
%     end
%     r(1)=sflux;
%     
%     % solve direct radiation
%     M1=diag(b)+diag(a,-1)+diag(c,1);
%     u=M1\r;
% 
% 
%     flux_u=mu0_local.*u(1:2:end);
%     flux_d=mu0_local.*u(2:2:end);
%     % save
%     flux_dn_direct1=flux_d;
%     flux_up_direct1=flux_u;
%     
%     %+++++++++++
%     % now diffuse
%     R=ref_diff;
%     T=trans_diff;
%     % build the matrix - a
%     for i=1:N
%         a(i*2-1)=1;
%         a(i*2)=-T(i);
%     end
%     a(end)=1;
%     % build the matrix - c
%     for i=1:N
%         c(i*2-1)=1;
%         c(i*2)=-T(i);
%     end
%     c(end)=1;
%     % build the matrix - b
%     for i=1:N
%         b(i*2)=-R(i);
%         b(i*2+1)=-R(i);
%     end
%     b(end)=-albedo;
% 
%     r(1)=0.;
%     r(end)=flux_dn_direct1(end).*albedo./mu0_local;
%     % N is the surface
%     for i=1:N-1
%         source_up(i)=flux_d(i)*ref_diff(i)+flux_u(i)*trans_dir_diff(i)/mu0_local;
%         source_dn(i)=flux_u(i)*ref_diff(i)+flux_d(i)*trans_dir_diff(i)/mu0_local;
%     end
%     source_up(N)
%     % build the matrix - r
%     for i=1:N
%         r(i*2)=source_up(i); %  source_up(i);
%         r(i*2+1)=source_dn(i); %source_dn(i);
%     end
%     
%     % solve diffuse radiation
%     M1=diag(b)+diag(a,-1)+diag(c,1);
%     u=M1\r;
% 
% 
%     flux_u=u(1:2:end);
%     flux_d=u(2:2:end);
% 
% end

albedo_at=zeros([N+1,1]);

if ~shortwave
    reflectance=R;
    transmittance=T;
    flux_dn=zeros([N+1,1]);
    flux_up=zeros([N+1,1]);
    source=zeros([N+1,1]);
    inv_denominator=zeros([N,1]);

    albedo_at(N+1)=albedo;
    
    % at the surface the source is thermal emission
    source(N+1)=emiss*pi*blt*tl(end)^4./(t(end)^4);
    
    % work bak up through the atmosphere and compute the albedo of the
    % entire earth / atmosphere system below that half-level and also the
    % source, which is the upwelling flux due to emission below that level
    for i=N:-1:1
        % lacis and hansen (1974) eq 33, shonk and hogan (2008) eq 10
        inv_denominator(i)=1/(1-albedo_at(i+1)*reflectance(i));
        % shonk and hogan (2008) eq 9, petty (2006) eq 13.81
        albedo_at(i)=reflectance(i)+...
            transmittance(i)*transmittance(i)*albedo_at(i+1)*...
            inv_denominator(i);
        % shonk and hogan (2008) eq 11
        source(i)=source_up(i)+transmittance(i)*(source(i+1)+...
            albedo_at(i+1)*source_dn(i))*inv_denominator(i);
    end
    % toa, no diffuse downwelling
    flux_dn(1)=0.;
    
    % toa, all upwelling is due to emission below that level
    flux_up(1)=source(1);
    
    % work back down through the atmosphere computing the fluxes at each
    % half-level
    for i=1:N
        % shonk and hogan (2008) eq 14 after simplification
        flux_dn(i+1)=(transmittance(i)*flux_dn(i)+reflectance(i)*source(i+1) ...
           +source_dn(i))*inv_denominator(i);
       
        % shonk and hogan (2008) eq 12
        flux_up(i+1)=albedo_at(i+1)*flux_dn(i+1)+source(i+1);
       
    end
    
end



if shortwave
    flux_dn_direct=zeros([N+1,1]);
    flux_dn_diffuse=zeros([N+1,1]);
    flux_up=zeros([N+1,1]);
    source=zeros([N+1,1]);
    inv_denominator=zeros([N,1]);
    % ref_diff          - reflectance
    % trans_diff        - transmittance
    % ref_dir           - ref_dir
    % trans_dir_diff    - trans_dif_diff
    % trans_dir_dir     - trans_dir_dir
    reflectance = ref_diff;
    transmittance = trans_diff;

    flux_dn_direct(1)=sflux;
    for i=1:N
        flux_dn_direct(i+1) = flux_dn_direct(i)*trans_dir_dir(i);
    end

    % diffuse source at the surface - direct beam reflected back into the
    % diffuse stream
    source(N+1)=albedo*flux_dn_direct(N+1)*cos_theta_s;

    % work back up through the atmosphere and compute the albedo of the eart /
    % atmosphere system below that half-level and also the source, which is the
    % upwelling flux due to the direct radiation that is scattered below that
    % level
    albedo_at(N+1)=albedo;
    for i=N:-1:1
        % lacis and hansen, eq 33; shonk and hogan (2008) eq 10
       inv_denominator(i) = 1/(1-albedo_at(i+1)*reflectance(i)); 
       % shonk and hogan (2008), eq 9, petty (2006) eq 13.81
       albedo_at(i)=reflectance(i)+...
           transmittance(i)*transmittance(i)*albedo_at(i+1)*inv_denominator(i);
       source(i)=ref_dir(i)*flux_dn_direct(i)+transmittance(i)*(source(i+1)+ ...
           albedo_at(i+1)*trans_dir_diff(i)*flux_dn_direct(i))*inv_denominator(i);
    end

    % at the toa thee is no diffuse downwelling radiation 
    flux_dn_diffuse(1)=0.;

    % at the toa, all upwelling radiation is due to scattering by the direct 
    % beam below that level
    flux_up(1) = source(1);

    % work back down through the atmosphere computing the fluxes at each
    % half-level
    for i=1:N
        % shonk and hogan (2008) eq 14 (after simplification)
        flux_dn_diffuse(i+1)=(transmittance(i)*flux_dn_diffuse(i)+...
           reflectance(i)*source(i+1)+trans_dir_diff(i)*flux_dn_direct(i)) ...
           * inv_denominator(i);
       % shonk and hogan (2008) eq 12
       flux_up(i+1)=albedo_at(i+1)*flux_dn_diffuse(i+1)+source(i+1);
       flux_dn_direct(i)=flux_dn_direct(i)*cos_theta_s;
    end
    flux_dn_direct(N+1)=flux_dn_direct(N+1)*cos_theta_s;
end






% flux_d(1:end)=flux_d(1:end)+cos_theta_s*sflux*exp(-[0; tauc(:); ]/cos_theta_s);

% if shortwave
%     % see Liou, 3.5.3 (and 3.5.1a), page 108
%     flux_dn_sw=flux_dn_direct+flux_dn_diffuse;
%     flux_up_sw=flux_up;
%     
%     plot(-diff(flipud(flux_up_sw(1:end))-flipud(flux_dn_sw(1:end)))./dz.*86400./1005,...
%     flipud(z),'linewidth',2)
% 
%     if false
%         tot=-diff(flipud(flux_up_sw(1:end))-flipud(flux_dn_sw(1:end)))./dz.*86400./1005 +...
%             -diff(flipud(flux_up)-flipud(flux_dn))./dz.*86400./1005;
%         plot(tot,flipud(z),'linewidth',2)
%     end
% else
%     % the IR heating rate - minus sign removed
%     % see Liou, 4.7.1, and 4.7.2, page 160
%     plot(-diff(flipud(flux_up)-flipud(flux_dn))./dz.*86400./1005,flipud(z),'linewidth',2)
% end

