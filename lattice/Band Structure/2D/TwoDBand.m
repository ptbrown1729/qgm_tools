function [ kvectsBZ,lambda ] = TwoDBand(potential,a1,a2)
%Takes potential (in joules) as argument.

%Solve matrix equation 
%Ec(k-K)=Sum_K'[{U(K'-K)+(hbar^2/2m)(k-K')^2*Delta(K,K')}c(k-K')]
%which we recast as
%Ec_k(K)=Sum_K'[{U(K,K')+(hbar^2/2m)*(k-K')^2*Delta(K,K')}c_k(K')]

%run('\\128.112.86.75\lithium\Simulation\Band Structure\init')

%fundamental constants
hbar=1.054e-34; %J.s
m=9.10938291e-31; %kg
c=2.99792458e8; %m/s
e=1.602176565e-19; %C

%Define lattice vectors in real space, in nm.
amax=max(norm(a1),norm(a2));
unit_cell_area=norm(cross([a1,0],[a2,0]));
%Set number of sites in each direction after which to impose periodic b.c.
N1=51;
N2=51;

%Define potential in Joules
potential_form=@(x,y) potential(x,y);

%Need a function which is our potential in the unit cell and zero
%outside the unit cell.
 xv=[0,a1(1),a1(1)+a2(1),a2(1)];
 yv=[0 a1(2), a1(2)+a2(2), a2(2)];
 charfn_unit_cell=@(x,y)inpolygon(x,y,xv,yv);
%potential_unit_cell=@(x,y)(potential_form(x,y).*WignerSeitz(x,y,a1,a2));

%Remove offset
offset=0;%integral2(@(x,y)(potential_form(x,y).*charfn_unit_cell(x,y)),-amax,amax,-amax,amax)/unit_cell_area;
potential=@(x,y)(potential_form(x,y)-offset).*charfn_unit_cell(x,y);

%Find reciprocal lattice basis vectors
[b1,b2]=getRecpBasis(a1,a2);
%Set number of Reciprocal lattice vectors to include in computation.
nrecb1=7;
nrecb2=7;

%Create Reciprocal lattice vectors
recb1_n=(-(nrecb1-1)/2:1:(nrecb1-1)/2);
recb1s=kron(transpose(recb1_n),b1);
recb2_n=(-(nrecb2-1)/2:1:(nrecb2-1)/2);
recb2s=kron(transpose(recb2_n),b2);

%create allowed k-vectors in first Brillouin zone.
kb1_n=(-(N1-1)/2:1:(N1-1)/2)/N1;
BZb1s=kron(transpose(kb1_n),b1);
kb2_n=(-(N2-1)/2:1:(N2-1)/2)/N2;
BZb2s=kron(transpose(kb2_n),b2);

%define complex exponential function
four_exp=@(x,y,k1,k2) exp(1i*(x*k1+y*k2));


%create potential matrix
%matrix indices=dictionary=(-n1,-n2),(-n1,-n2+1),...(-n1,n2),(-n1+1,-n2)...
%this means that the index(l1,l2)=(n1+l1)*(2n2+1)+(l2+n2+1)
%implement this order using kronecker product as below.


recpvects=kron(recb1s,ones(1,length(recb2_n)).')+kron(ones(1,length(recb1_n)).',recb2s);
Umat=zeros(nrecb1*nrecb2);
for ii=1:size(recpvects,1)
    for jj=1:size(recpvects,1)
        Q=recpvects(jj,:)-recpvects(ii,:);
        Umat(ii,jj)=(1/unit_cell_area)*integral2(@(x,y)(potential(x,y).*four_exp(x,y,Q(1),Q(2))),-amax,amax,-amax,amax);
    end
end

%create array to hold eigenvalues.
lambda=[];
kvects=kron(ones(N2,1),BZb2s)+kron(BZb1s,ones(N1,1));
for ii=1:size(kvects,1)
    k=kvects(ii,:);
    delta_vect=zeros(1,nrecb1*nrecb2);
    for jj=1:size(recpvects,1)
        Q=recpvects(jj,:);
        delta_vect(jj)=(hbar^2/(2*m)).*norm(k-Q)^2;
    end
    delta_part=diag(delta_vect);
    %Create and diagonalize matrix. Append as next row of lambda
    FullMatrix=Umat+delta_part;
    lambda=cat(1,lambda,transpose(eig(FullMatrix)));
end

%reduce k's to Brillouin zone
kvalscell=num2cell(kvects,2);
getkBZ=@(pt)reduceToWignerSeitz(pt(1),pt(2),b1,b2);
kvectsBZ=cell2mat(cellfun(getkBZ,kvalscell,'UniformOutput',0));

end

