function [qvects,Es,psiKs,recpvects] = OneDBand(potential,a)
%[qvects,Es,psiKs,recpvects] = OneDBand(potential,a)
%Use SI units
%qvectors is a 1xNbrill matrix
%bvects is a 1xNrec matrix
%lambda, a single column is a band.
%psiKs. Column is a wavefunction. 2D matrix of wavefunctions for given q.
%Third dimension changing q.

%For a 1D Hamiltonian we can find a real solution. However, these
%eigenstates may not be eigenstates of the translation operator.

%PsiKs associated with different q's are automatically orthogonal due to
%the fact no q is commensurate with a lattice vector. Even though they are
%all returned as vectors of the same lengths, these live in different
%spaces for different q's.

%Solve matrix equation 
%E*c(k-K)=Sum_K'[{U(K'-K)+(hbar^2/2m)(k-K')^2*Delta(K,K')}*c(k-K')]
%which we recast as
%E*c_k(K)=Sum_K'[{U(K,K')+(hbar^2/2m)*(k-K')^2*Delta(K,K')}*c_k(K')]

run('constants.m')

potential_form=@(r)potential(r);
%((heaviside(-r+a/2).*(20*r/a)+heaviside(r-a/2).*(20-20*r/a)))*e;
%want our potential to be zero on average.
%offset=integral(@(r)potential_form(r.*a),0,1);
offset=integral(@(r)potential_form(r.*a),-0.5,0.5);
potential=@(r)potential_form(r)-offset;

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

%define the complex exponential function we'll need for our fourier
%integrals
four_exp=@(r,k) exp(1i*k*r);


Umat=zeros(Nrec);
for ii=1:Nrec
    for jj=1:Nrec
        Q=recpvects(jj)-recpvects(ii);
        %Umat(ii,jj)=integral(@(r)(potential(r.*a).*four_exp(r.*a,Q)),0,1);
        Umat(ii,jj)=integral(@(r)(potential(r.*a).*four_exp(r.*a,Q)),-0.5,0.5);
        %choose the symmetric BZ, since then get real coefficients.
    end
end

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
end

