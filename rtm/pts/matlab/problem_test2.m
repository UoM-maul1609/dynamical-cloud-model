function [u1,u2,u3]=problem_test2
% serial algorithm 1
sz=100;
A=2.*eye(sz,sz)+... % diagonal
    [zeros(sz,1),[eye(sz-1,sz-1);zeros(1,sz-1)]]+... % upper
    [[zeros(1,sz-1);0.5.*eye(sz-1,sz-1)],zeros(sz,1)];    % lower
b=ones(sz,1);

u1=A\b;

% serial algorithm 2
b=2.*ones(sz,1);
a=0.5.*ones(sz-1,1); % lower diagonal
a(50)=0.2;
c=ones(sz-1,1); % upper diagonal
r=ones(sz,1); % result

u2=tridag_ser(b,a,c,r,sz);


% parallel algorithm https://web.alcf.anl.gov/~zippy/publications/partrid/partrid.html
p=4;
assert(sz./p==round(sz./p));
u3=tridag_par(b,a,c,r,sz,p);


function u=tridag_ser(b,a,c,r,sz)

bet=b(1);
u=zeros(sz,1);
gam=zeros(sz,1);
u(1)=r(1)./bet;
if(bet==0) return;end
for j=2:sz
    gam(j)=c(j-1)./bet;
    bet=b(j)-a(j-1).*gam(j);
    if(bet==0) return;end;
    u(j)=(r(j)-a(j-1).*u(j-1))./bet;
end

for j=sz-1:-1:1
    u(j)=u(j)-gam(j+1).*u(j+1);
end

function u=tridag_par(bn,an,cn,rn,sz,p)

% make sub-arrays+++++++++++++++++++++++++
m=sz./p;
a=zeros(m,p);
b=zeros(m,p);
c=zeros(m,p);
r=zeros(m,p);
for i=1:p
    b(:,i)=bn(1+(i-1).*m:i.*m);
    if(i==1) % lower diagonal
        a(1,i)=0;
        a(2:m,i)=an(1:m-1);
    else
        a(:,i)=an([(i-1).*m:i.*m-1]);
    end
    
    if(i==p) % upper diagonal
        c(1:m-1,i)=cn([1+(i-1)*m:i*m-1]);
        c(m,i)=0.;
    else
        c(:,i)=cn(1+(i-1).*m:i.*m);        
    end
    
    r(:,i)=rn(1+(i-1).*m:i.*m);
end
%------------------------------------------



% part 1: forward elim, back sub and forward sub
xuh=zeros(m,p);
xlh=zeros(m,p);
xr=zeros(m,p);
for j=1:p
    % forward elimination
    xuh(1,j)=c(1,j)/b(1,j);
    xlh(1,j)=r(1,j)/b(1,j);
    for i=2:m
        denom=b(i,j)-a(i,j)*xuh(i-1,j); % mistake had a as c here
        assert(denom~=0);
        xuh(i,j)=c(i,j)/denom;
        xlh(i,j)=(r(i,j)-a(i,j)*xlh(i-1,j))/denom;
    end
    % back substitution
    xr(m,j)=xlh(m,j);
    xlh(m,j)=-xuh(m,j);
    xuh(m,j)=a(m,j)/b(m,j);
    for i=m-1:-1:1
        xr(i,j)=xlh(i,j)-xuh(i,j)*xr(i+1,j);
        xlh(i,j)=-xuh(i,j)*xlh(i+1,j);
        denom=b(i,j)-c(i,j)*xuh(i+1,j);
        assert(denom ~=0);
        xuh(i,j)=a(i,j)/denom; % minus sign wrong in paper
    end
    
    % forward substitution
    xuh(1,j)=-xuh(1,j);
    for i=2:m
        xuh(i,j)=-xuh(i,j)*xuh(i-1,j);
    end
    
end


% part 2: send information back
log2P=log(p)./log(2);
OutData=zeros(8*2^log2P,p);
OutData(1,:)=-1;
OutData(2,:)=xuh(1,:);
OutData(3,:)=xlh(1,:);
OutData(4,:)=-xr(1,:);
OutData(5,:)=xuh(m,:);
OutData(6,:)=xlh(m,:);
OutData(7,:)=-1;
OutData(8,:)=-xr(m,:);

% OutDataLocal=reshape(OutData,[8*p 1]);
% OutDataLocal=zeros(8*2^log2P,1);
% OutData
% for j=1:p
%     for i=0:log2P-1
%         nxfer=8*(2^i);
%         ToProc=1 + mod(j+1-2^i+2*p,p);
%         FromProc=1 + mod(j+2^i,p);
%         [ToProc, FromProc]
%         OutData([1:8]+(FromProc-1)*8,ToProc)=OutData(1:8,FromProc);
%     end
% end
% OutData
d=OutData(1:8,:);

OutData=repmat(d(:),[1 p]);
if(p==1)
    u=xr;
    return;
end

% put outdata into reduced tridiagonal form:
nsig=8*p;
ifirst=8*(p-[p:-1:1])+5
sz2=2*p-2;
x=zeros(m,p);
for j=1:p
    for i=1:2*p-2
        ibase=mod(ifirst(1)+4*(i-1),nsig)
        reduca(i,j)=OutData(ibase,1);
        reducb(i,j)=OutData(ibase+1,1);
        reducc(i,j)=OutData(ibase+2,1);
        reducr(i,j)=OutData(ibase+3,1);
    end
    % solve reduced system
    coeffs(:,j)=tridag_ser(reducb(:,j),reduca(2:end,j),...
        reducc(1:end-1,j),reducr(:,j),sz2);
end


% pick out the appropriate elements of coeffs
for j=1:p
    if(j ~= 1)
        uhcoeff=coeffs(2*j-2,j);
    else
        uhcoeff=0;
    end
    
    if(j ~= p)
        lhcoeff=coeffs(2*j-1,j);
    else
        lhcoeff=0;
    end
    % compute the final solution
    for i=1:m
        x(i,j)=xr(i,j)+uhcoeff*xuh(i,j)+lhcoeff*xlh(i,j);
    end
end



u=reshape(x,[m*p 1]);





