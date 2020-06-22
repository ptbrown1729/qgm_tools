%[qvects,Es,psiKs,recpvects] = OneDBand(potential,a)
%Use SI units
%qvectors is a 1xNbrill matrix
%bvects is a 1xNrec matrix
%lambda, a single column is a band.
%psiKs. Column is a wavefunction. 2D matrix of wavefunctions for given q.
%Third dimension changing q.

%PsiKs associated with different q's are automatically orthogonal due to
%the fact no q is commensurate with a lattice vector. Even though they are
%all returned as vectors of the same lengths, these live in different
%spaces for different q's.

%Solve matrix equation 
%E*c(k-K)=Sum_K'[{U(K'-K)+(hbar^2/2m)(k-K')^2*Delta(K,K')}*c(k-K')]
%which we recast as
%E*c_k(K)=Sum_K'[{U(K,K')+(hbar^2/2m)*(k-K')^2*Delta(K,K')}*c_k(K')]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set klatt and depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Potential = Depth/2*cos(kx)
%Er = hbar^2*klight^2/(2*mLi);

klatt = 2*pi/1064e-9*sqrt(2);
Er = hbar^2 / (2*mLi) * (klatt / 2)^2;
Depth = 6*Er;
a = 2*pi/klatt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define useful constants
 c = 299792458;
 mLi = 9.9883e-27;
 h = 6.6262e-34;
 hbar = h/(2*pi);

%Set number of sites after which impose periodic b.c.'s
%Ensure this is odd.
Nrec=201; %dimension of matrices
Nbrill=201; %number of matrices

%Reciprocal vectors to compute FT with
b=2*pi/a;
%Reciprical vectors for matrix indices.
recpvects=transpose(-(Nrec-1)/2*b:b:(Nrec-1)/2*b);
%Create allowed k vectors in the first Brillouin zone.
qvects=(2*pi/(Nbrill*a))*(-(Nbrill-1)/2:1:(Nbrill-1)/2);


%%%Uncomment this part for a more general potential. For simple enough
%%%potentials it is probably better to do fourier transform by hand
%potential_form=@(r)potential(r);
%want our potential to be zero on average.
%offset=integral(@(r)potential_form(r.*a),0,1);
%offset=integral(@(r)potential_form(r.*a),-0.5,0.5);
%potential=@(r)potential_form(r)-offset;
%define the complex exponential function we'll need for our fourier
%integrals
%four_exp=@(r,k) exp(1i*k*r);
% Umat=zeros(Nrec);
% for ii=1:Nrec
%     for jj=1:Nrec
%         Q=recpvects(jj)-recpvects(ii);
%         %Umat(ii,jj)=integral(@(r)(potential(r.*a).*four_exp(r.*a,Q)),0,1);
%         Umat(ii,jj)=integral(@(r)(potential(r.*a).*four_exp(r.*a,Q)),-0.5,0.5);
%         %choose the symmetric BZ, since then get real coefficients.
%     end
% end

%For potential of the form 2*Vo*cos(kx)
Umat = Depth/4*(diag(ones(1,Nrec-1),1)+diag(ones(1,Nrec-1),-1));


%Create empty array to hold eigenvalues. Each column is a single band. 
Es=[];
psiKs=[];
%sum over all k-vectors which are allowed.
for  ii=1:Nbrill
    q=qvects(ii);
    %get diagonal parts of matrix
    delta_vect=zeros(1,Nrec);
    for jj=1:Nrec
        Q=recpvects(jj);
        delta_vect(jj)=(hbar^2/(2.*mLi)).*norm(q-Q)^2;
    end
    delta_part=diag(delta_vect);
    %create and diagonalize matrix
    FullMatrix=Umat+delta_part;
    %append as the next row of lambda
    [psiK_current,lambda_current] = eig(FullMatrix);
    lambda_current=diag(lambda_current);
    Es=cat(1,Es,transpose(lambda_current));
    %lambda=cat(1,lambda,transpose(eig(FullMatrix)));
    psiKs=cat(3,psiKs,psiK_current);
   
    
end

%for inversion symmetric potential, we can also pick a smooth gauge. We can
%do this in K-space as well as real space.

%pick some element, and compare the signs of it.
Kindex = floor(size(psiKs,1)/2);
Ksigns = -1*real(sign(psiKs(Kindex,:,:)));
valsMat = repmat(Ksigns,[Nrec,1,1]);
psiKs = psiKs.*valsMat;

figname = sprintf('Band for %0.1f Er deep lattice',Depth/Er);
figure('name',figname);
NBandsToPlot = 3;
for ii = 1:NBandsToPlot
    plot(Es(:,ii)/Er)
    hold on;
end

ylabel('Energy (Er)')

fprintf('Depth = %0.1fEr \n',Depth/Er);
fprintf('t = %0.5f Er \n', 0.25 * (max(Es(:,1))-min(Es(:,1))) / Er);
fprintf('Band Gap 1 to 2 = %0.1f KHz \n',(min(Es(:,2))-max(Es(:,1)))/h/1e3);
fprintf('Band Gap 2 to 3 = %0.1f KHz \n',(min(Es(:,3))-max(Es(:,2)))/h/1e3);
fprintf('Band Gap 1 to 3 = %0.1f KHz \n',(min(Es(:,3))-max(Es(:,1)))/h/1e3);

