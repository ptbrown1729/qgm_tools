function qbz = reduceToWignerSeitz(qx,qy,b1,b2)
%reduceToBZ
%takes a point and two reciprocal lattice vectors
%and finds an equivalent point in the first Brillouin zone.


testBZ=@(pt)WignerSeitz(pt(1),pt(2),b1,b2);


%strategy b1 and b2 are independent, so we can find real number a*b1+c*b2=q.
M=[b1;b2].';
%Shouldn't use inv...but having trouble with \
q=[qx,qy];
soln=inv(M)*q.';
%now we want to find the closest integer values. There are only four
%choices remaining.
n1=floor(soln(1));
n2=ceil(soln(1));
m1=floor(soln(2));
m2=ceil(soln(2));

qbz=0;
for n=[n1 n2]
    for m=[m1 m2]
        if testBZ(q-n*b1-m*b2)==1
          qbz=q-n*b1-m*b2;
        else
        end
    end
end


end

