%Get tunneling values using real space diagonalization of lattice
%Hamiltonian plus wannier1D.

run('constants.m')

%define constants
lambda = 1064e-9;
k=2*pi/lambda;
Er=hbar^2*k^2/(2*mLi);
V1=@(x)0.2*Er*cos(2*k*x);
V2=@(x)0.7*Er*cos(2*k*x);
V3=@(x)1.2*Er*cos(2*k*x);
V4=@(x)3*Er*cos(2*k*x);
V5=@(x)5*Er*cos(2*k*x);
V6=@(x)10*Er*cos(2*k*x);

[H1,xpts] = Hlattice(V1,lambda/2);
H2 = Hlattice(V2,lambda/2);
H3 = Hlattice(V3,lambda/2);
H4 = Hlattice(V4,lambda/2);
H5 = Hlattice(V5,lambda/2);
H6 = Hlattice(V6,lambda/2);

wanns = zeros(length(H1),51,6);
wanns(:,:,1) = Wannier1D(H1,xpts,51);
wanns(:,:,2) = Wannier1D(H2,xpts,51);
wanns(:,:,3) = Wannier1D(H3,xpts,51);
wanns(:,:,4) = Wannier1D(H4,xpts,51);
wanns(:,:,5) = Wannier1D(H5,xpts,51);
wanns(:,:,6) = Wannier1D(H6,xpts,51);


tnn = [wanns(:,20,1)'*H1*wanns(:,21,1),wanns(:,20,2)'*H2*wanns(:,21,2),wanns(:,20,3)'*H3*wanns(:,21,3),wanns(:,20,4)'*H4*wanns(:,21,4),wanns(:,20,5)'*H5*wanns(:,21,5),wanns(:,20,6)'*H6*wanns(:,21,6)];
tnnn = [wanns(:,20,1)'*H1*wanns(:,22,1),wanns(:,20,2)'*H2*wanns(:,22,2),wanns(:,20,3)'*H3*wanns(:,22,3),wanns(:,20,4)'*H4*wanns(:,22,4),wanns(:,20,5)'*H5*wanns(:,22,5),wanns(:,20,6)'*H6*wanns(:,22,6)];
tnnnn = [wanns(:,20,1)'*H1*wanns(:,23,1),wanns(:,20,2)'*H2*wanns(:,23,2),wanns(:,20,3)'*H3*wanns(:,23,3),wanns(:,20,4)'*H4*wanns(:,23,4),wanns(:,20,5)'*H5*wanns(:,23,5),wanns(:,20,6)'*H6*wanns(:,23,6)];
figure
semilogy([0.2,0.7,1.2,3,5,10],abs(tnn/Er),'r.');
hold on;
semilogy([0.2,0.7,1.2,3,5,10],abs(tnnn/Er),'b.');
hold on;
semilogy([0.2,0.7,1.2,3,5,10],abs(tnnnn/Er),'g.');
legend('t nn','t nnn','t nnnn');
title('Tunneling vs. Lattice Depth')
xlabel('Lattice Depth [Er]')
ylabel('Tunneling [Er]')

figure
colorset=jet(6);
for ii=1:6
plot(xpts,wanns(:,20,ii),'color',colorset(ii,:))
hold on;
end
legend('0.2','0.7','1.2','3','5','10');