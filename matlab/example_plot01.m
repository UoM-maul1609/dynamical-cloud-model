% example plotting script
print_image=false;
nc=netcdf('/tmp/output.nc');

imax=length(nc{'time'}(:));
n=1;
for i=1:imax
    pcolor(nc{'y'}(:)./1000,nc{'z'}(:)./1000,...
        squeeze(nc{'th'}(i,1,:,:))');shading flat
    caxis([-0.1 0.1]);
    xlabel('y (km)');ylabel('z (km)');
    h=colorbar;
    ylabel(h,'\theta (K)');
    
    pause(0.1);
    if print_image
        if(n==1)
            mkdir /tmp/pics/
        end
        eval(['print -dpng /tmp/pics/dcm_output_', ...
            num2str(n,'%03d'),'.png']);
    end
    n=n+1;
end

close(nc);
