%fundamental constants
hbar=1.054e-34; %J.s
m=9.10938291e-31; %kg
c=2.99792458e8; %m/s
e=1.602176565e-19; %C

%Set time periodicity
omega=2*pi*6e16;
T=2*pi/omega; 
%set lattice periodicity
a=0.2e-9; %in m
%Define potential, in J.
k_lattice=2*pi/a;
Erecoil=hbar^2*k_lattice^2/(2*m);
V=15;
b=0.0/k_lattice;
potential_form=@(x,t)V*cos(k_lattice*x+b.*cos(omega*t)).^2*e;

%select number of Fourier modes, Recp Vects, and Brillouin zone q vects to use.
%These all should be odd.
Nf=3; 
Nrec=31; %Nf*nrec is dimension of matrices
Nbrill=301; %number of matrices to diagonalize

%create these
nfvals=-(Nf-1)/2:1:(Nf-1)/2;
RecBaseVect=2*pi/a;
RecVects=-(Nrec-1)/2*RecBaseVect:RecBaseVect:(Nrec-1)/2*RecBaseVect;
Qvectors=(2*pi/(Nbrill*a))*(-(Nbrill-1)/2:1:(Nbrill-1)/2);

%define the complex exponential function we'll need for our fourier
%integrals
four_space_exp=@(x,k) exp(1i*x*k);
four_time_exp=@(t,w) exp(1i*w*t);

%create matrices
Umat=zeros(Nrec,Nrec,Nf);
for nn=1:Nf
    for ii=1:Nrec
        for jj=1:Nrec
            Q=RecVects(jj)-RecVects(ii);
            angfreq=omega*nfvals(nn);
            Umat(ii,jj,nn)=integral2(@(u,v)(potential_form(u.*a,v.*T).*four_space_exp(u,Q.*a).*four_time_exp(v,angfreq.*T)),0,1,0,1);
        end
    end
end

%Create skeleton matrix. This contains terms that do not vary with q. 
H_skeleton=kron(diag(nfvals),eye(Nrec))*hbar*omega;
%Blocks for each time fourier component
for nn=1:Nf
    H_skeleton=H_skeleton+kron(diag(ones(1,Nf-abs(nfvals(nn))),nfvals(nn)),Umat(:,:,nn));
end

%create array to hold eigenvalues for each q value
lambda=[];
for ii=1:Nbrill
    k=Qvectors(ii);
    %delta part
    delta_vect=zeros(1,Nrec);
    for jj=1:Nrec
        Q=RecVects(jj);
        delta_vect(jj)=(hbar^2/(2*m)).*norm(k-Q)^2;
    end
    delta_part=kron(eye(Nf),diag(delta_vect));
    %create and diagonalize matrix. Append as next row of lambda
    floquet_H=delta_part+H_skeleton;
    lambda=cat(1,lambda,transpose(eig(floquet_H)));   
end

%plot the band structure over the first Brillouin zone.
subplot(2,2,1)
%plot(Kvectors/(pi/a),lambda)
plot(Qvectors/(pi/a),lambda(:,1:9)/e)
xlabel('Crystal Momentum [pi/a]')
ylabel('Energy [eV]')
subplot(2,2,2)
x=linspace(0,1,100);
plot(x,potential_form(a*x,0)/e,'.')
xlabel('Lattice Spacing [a]')
ylabel('Potential [eV]')
    