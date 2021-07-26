lambdas=linspace(100e-9,100e-6,10000);
average=true;
% solar_cutoff=5e-6;
if average
    for i=1:length(lambdas)-1
        [blt(i),x_upper(i),x_lower(i),sum_upper(i),sum_lower(i)]=...
            inline_planck(lambdas(i),lambdas(i+1),290);
    end
end

result=zeros(length(lambdas),1);
result2=zeros(length(lambdas),1);
qc=0.3e-3;

mu=7;
Dm=20e-6;
% lam=(pi/6*1000*gamma(mu+4)/gamma(mu+1)*nc/qc)^(1/3);
lam=gamma(mu+2)/gamma(mu+1)/Dm;



n0=qc*lam^(mu+4)/gamma(mu+4)/(pi./6*1000);

refractive_indices;
% lambdas=10e-8;
% lam=603223.31540111464;
% n0=1.0970550231289439E+052;
% mu= 7.3163536020529483;
for i=1:length(lambdas)
    mtun=0.5;
    lambda= lambdas(i); 
    if(lambda<min(lam_h2o)) 
        nrwbin=nr_h2o(1); 
        niwbin=ni_h2o(1); 
    elseif(lambda>max(lam_h2o))
        nrwbin=nr_h2o(end); 
        niwbin=ni_h2o(end); 
    else
        nrwbin=interp1(lam_h2o,nr_h2o,lambda);
        niwbin=interp1(lam_h2o,ni_h2o,lambda);
    end
    sigma=pi/4;
    delta=2;
    g=8*pi*niwbin/(3*lambda);
    e0=0.25+0.6*(1.-exp(-8.*pi*niwbin/3.))^2; % eq 8
    e=e0/lambda; % equation 28
    ra=0.7393*nrwbin-0.6069; % equation 6
    rext=ra*0.5; % equation 11
    a1=0.25+0.25*exp(-1167.*niwbin); % equation5
    kmax=mtun/e0; % equation 9
    a2=ra/(kmax^mtun*exp(-mtun)*lambda^mtun);
    a3=rext/(kmax^mtun*exp(-mtun)*lambda^mtun);
    a4=0.06*pi/lambda;
    a5=(pi/lambda)^(-2./3.);
    a6=1.;

    gam_a_1=gamma(mu+3);
    gam_a_2=gamma(3+mtun+mu);

    babs= ...
        sigma*n0 * ...
            gam_a_1 / ...
            (lam^(3.+mu)) - ...
        sigma*n0 * ...
            gam_a_1 / ...
            ((lam+g)^(3.+mu)) + ...
        a1*sigma*n0 * ...
            gam_a_1 / ...
            ((lam+g)^(3.+mu)) - ...
        a1*sigma*n0 * ...
            gam_a_1 / ...
            ((lam+2.*g)^(3.+mu)) + ...
        a2*sigma*n0 * ...
            gam_a_2 / ...
            ((lam+e)^(3.+mu+mtun)) - ...
        a2*sigma*n0 * ...
            gam_a_2 / ...
            ((lam+e+g)^(3.+mu+mtun));
    m2=sigma*n0*gamma(mu+3)/lam^(mu+3);
    gam_3=gamma(mu+7/3);
    gam_4=gamma(mu+2);
    gam_5=gamma(mu+1);
    gam_6=gamma(mtun+mu+2);
    gam_7=gamma(mtun+mu+1);
    imag1=complex(0,1);
    ncom=complex(nrwbin,-niwbin);
    q=imag1*(ncom-complex(1,0))*complex(2*pi/lambda,0);
    
    % extinction
    test= ...
        pi*n0*gam_a_1 / (2.*lam^(3.+mu)) + ... % checked
        a3*pi*n0*gam_a_2 / (2.*(lam+e)^(3.+mu+mtun)) + ... % checked
        a6*sigma*a5*n0*gam_3*...
            (lam^(-(mu+7./3.))-(lam+a4)^(-(mu+7./3.)))+... % checked
        pi*n0*real( ...
         complex(gam_4,0.) / (q*(complex(lam,0.)+q)^(complex(mu+2.,0.))) +...
         complex(gam_5,0.) / (q*q)*((complex(lam,0.)+q)^(-(complex(mu+1.,0.))) - ...
            complex(lam,0.)^(-complex(mu+1.,0.))))+...
        a3*pi*n0*real( ...
         complex(gam_6,0.) / (q*(complex(lam+e,0.)+q)^ complex(mu+mtun+2.,0.)) +...
         complex(gam_7,0.) / (q*q)*((complex(lam+e,0.)+q)^(-complex(mu+mtun+1.,0.)) - ...
           complex(lam+e,0.)^(-complex(mu+mtun+1.,0.))));   

    result(i)=babs/m2;
    result2(i)=test/m2;
end

if average
    sigma_a_clouds=sum(result(1:end-1).*blt'.*...
        (lambdas(2:end)-lambdas(1:end-1))')./ ...
        sum(blt'.*(lambdas(2:end)-lambdas(1:end-1))');
    sigma_s_clouds=sum(result2(1:end-1).*blt'.*...
        (lambdas(2:end)-lambdas(1:end-1))')./ ...
        sum(blt'.*(lambdas(2:end)-lambdas(1:end-1))');
end
