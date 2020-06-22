run('constants.m')
Er = hbar^2*(2*pi/1064e-9)^2/(2*mLi)/2;
Depths = [7]; %in units of Er
Neigs = 1;
Nplots = Neigs + 2;
Ncols = ceil(sqrt(Nplots));
Nrows = ceil((Nplots)/Ncols);

% for ii = 1:length(Depths)
%     [H,V,xpts] = HlatticeWithTrap(2*pi*100,Depths(ii)*Er,532e-9);    
%     [psivects,Es] = getEigs1D(H,Neigs);
%     
%     figname = sprintf('Depth = %0.2f Er',Depths(ii));
%         figure('name',figname)
%     for jj = 1:Neigs
%         subplot(Nrows,Ncols,jj)
%         plot(psivects(:,jj));
%     end
%     subplot(Nrows,Ncols,Nplots-1)
%     plot(xpts,V(xpts))
%     subplot(Nrows,Ncols,Nplots)
%     plot(xpts,sum(psivects.^2,2))
% end

for ii = 1 : length(Depths)
     [H, V, xpts] = HlatticeWithTrap( 2*pi*100, Depths(ii) * Er, 532e-9 * sqrt(2) );    
     [psivects, Es] = getEigs1D(H, Neigs);
     
    
     
     figure()
     
     % plot density
     subplot(2, 1, 1);
     plot(xpts, sum(psivects.^2, 2) )
     title( sprintf('Depth = %0.2f Er', Depths(ii)) );
     
     % plot potential
     subplot(2, 1, 2)
     plot(xpts, V(xpts));
     xlabel('position');
     ylabel('V(x)');
     title('potential');
end

