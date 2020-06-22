function [H,V,x1pts,x2pts] = H_TwoParticles_DeltaFunction_HarmonicTrap(P)
% = Hgaussian2D(P)
%P = [a,Omega]
%H, the position space Hamiltonian is returned using a 'dictionary' basis.
%i.e. an eigenvector is made up of basis vectors in the following order...
%(x1,y1),(x1,y2),...(x1,yn),(x2,y1),...(xn,yn)
a = P(1);
Omega = P(2);

%useful constants
c = 299792458;
abohr = 5.2918e-11;
hbar = 1.0546e-34;
mLi = 9.9883e-27;
achar=sqrt(hbar/(2*mLi*Omega));

%Potential
%g = 4*pi*hbar^2*a/(mLi)*(1/(sqrt(2*pi)*a))^2; %units have to match 1D hamiltonian...this corresponds to integrating out other coords...
g = 4*pi*hbar^2*a/(mLi);
W = achar/5; %width of 'delta function'
V = @(x1,x2) g*(heaviside((x1-x2)+W/2).*heaviside((x2-x1))+W/2)/W + 0.5*mLi*Omega^2*(x1.^2+x2.^2); %some approximation of a delta function...

%create x-grid
EndPtsNum = 10;
Division = 20;
x1start = -EndPtsNum*achar;
x1end = EndPtsNum*achar;
dx1 = achar/Division;
x1pts = x1start:dx1:x1end;
Nx1 = length(x1pts);
%create y-grid
x2start = -EndPtsNum*achar;
x2end = EndPtsNum*achar;
dx2 = achar/Division;
x2pts = x2start:dx2:x2end;
Nx2 = length(x2pts);


%build position vector
pos = [kron(x1pts',ones(length(x2pts),1)),kron(ones(length(x1pts),1),x2pts')];
%(x1,y1),(x1,y2),...(x1,yn),(x2,y1),...
%Build Hamiltonian

Vmatr = sparse(1:Nx1*Nx2,1:Nx1*Nx2,V(pos(:,1),pos(:,2)));

Dx1Form = -2*sparse(1:Nx1,1:Nx1,ones(1,Nx1))+sparse(1:(Nx1-1),2:Nx1,ones(1,Nx1-1),Nx1,Nx1)+sparse(2:Nx1,1:(Nx1-1),ones(1,Nx1-1),Nx1,Nx1);
Dx1Sqrd = kron(Dx1Form,sparse(1:Nx2,1:Nx2,ones(1,Nx2)))/dx1^2;
Dx2Form = -2*sparse(1:Nx2,1:Nx2,ones(1,Nx2))+sparse(1:(Nx2-1),2:Nx2,ones(1,Nx2-1),Nx2,Nx2)+sparse(2:Nx2,1:(Nx2-1),ones(1,Nx2-1),Nx2,Nx2);
Dx2Sqrd= kron(sparse(1:Nx1,1:Nx1,ones(1,Nx1)),Dx2Form)/dx2^2;

H = -0.5*hbar^2/(2*mLi)*(Dx1Sqrd+Dx2Sqrd) + Vmatr;




end

