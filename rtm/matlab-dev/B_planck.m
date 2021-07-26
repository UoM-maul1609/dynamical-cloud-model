function B_planck

lambda=[0.105:0.005:1].*1e-6;
h=6.6256e-34;
c=2.9979e8;
kb=1.38e-23;
T=5778;
B=2.*h.*c.^2./(lambda.^5.*(exp(h.*c./(lambda.*kb.*T))-1));

figure;
plot(lambda,B);
% total=5.67e-8.*T.^4;

% f=inline('pi.*1370*2.*6.6256e-34.*2.9979e8.^2./(lambda.^5.*(exp(6.6256e-34.*2.9979e8./(lambda.*1.38e-23.*5778))-1))./(5.67e-8*5778.^4)','lambda');
f=inline('2.*6.6256e-34.*2.9979e8.^2./(lambda.^5.*(exp(6.6256e-34.*2.9979e8./(lambda.*1.38e-23.*5778))-1))','lambda');
total=quad(f,1e-9,10e-6,1e-8);
B=B./total.*1370;
Btoa=B'.*1e-6;

for i=1:length(Btoa)
    fprintf('%.3f %4f\n',lambda(i).*1e6,Btoa(i));
end


