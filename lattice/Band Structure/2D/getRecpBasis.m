function [b1,b2] = getRecpBasis(a1,a2)
%Generate the reciprocal basis vectors from the real space
%basis vectors.

%write the conditions bi.aj=2pi delta_ij in matrix form A[b1 b2]=[2pi 0 0 2pi]
A=zeros(4,4);
A(1,1)=a1(1);
A(1,2)=a1(2);
A(2,3)=a1(1);
A(2,4)=a1(2);
A(3,1)=a2(1);
A(3,2)=a2(2);
A(4,3)=a2(1);
A(4,4)=a2(2);
%solution vector
vect=[2*pi 0 0 2*pi];
%invert matrix. 
b=vect/A;
b1=b(1:2);
b2=b(3:4);



end

