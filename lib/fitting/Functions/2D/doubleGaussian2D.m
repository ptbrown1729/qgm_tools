function [V,PNames,ArgStr] = doubleGaussian2D(P,X,Y)
%P = [CxG1,CyG1,SxG1,SyG1,AmpG1,ThetaG1,CxG2,CyG2,SxG2,SyG2,AmpG2,ThetaG2,Bg]
%Rotated 2D Gaussians
PNames = {'CxG1','CyG1','SxG1','SyG1','AG1','ThetaG1','CxG2','CyG2','SxG2','SyG2','AG2','ThetaG2','Bg'};
ArgStr = '';

if isempty(P) || isempty(X) || isempty(Y)
    V = 0;
else
%         CxG1 = P(1); CyG1 = P(2); SxG1 = P(3); SyG1 = P(4); AG1 = P(5); ThetaG1 = P(6);
%         CxG2 = P(7); CyG2 = P(8); SxG2 = P(9); SyG2 = P(10); AG2 = P(11); ThetaG2 = P(12);
%         Bg = P(13);
%         V = Bg...
%             +abs(AG1)*(exp(-0.5*((X-CxG1)*cos(ThetaG1)-(Y-CyG1)*sin(ThetaG1)).^2./(SxG1^2)...
%             -0.5*((Y-CyG1)*cos(ThetaG1)+(X-CxG1)*sin(ThetaG1)).^2./(SyG1^2)))...
%             +abs(AG2)*(exp(-0.5*((X-CxG2)*cos(ThetaG2)-(Y-CyG2)*sin(ThetaG2)).^2./(SxG2^2)...
%             -0.5*((Y-CyG2)*cos(ThetaG2)+(X-CxG2)*sin(ThetaG2)).^2./(SyG2^2)));
        V = gaussian2D([P(1:6), 0], X, Y) + gaussian2D(P(7:13), X, Y);

end
end

