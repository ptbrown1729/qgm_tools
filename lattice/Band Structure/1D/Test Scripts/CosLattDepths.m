%Get band structures for different lattice depths. Use k-space functions.
lambda=1064e-9;
k=2*pi/lambda;
Er=hbar^2*k^2/(2*mLi);
V1=@(x)0.2*Er*cos(2*k*x);
V2=@(x)0.7*Er*cos(2*k*x);
V3=@(x)1.2*Er*cos(2*k*x);
V4=@(x)3*Er*cos(2*k*x);
V5=@(x)5*Er*cos(2*k*x);
V6=@(x)10*Er*cos(2*k*x);

[qs,Es]=OneDBand(V1,lambda/2);
subplot(2,3,1)
plot(qs,Es(:,1)/Er)
hold on;
plot(qs,Es(:,2)/Er)
hold on;
plot(qs,Es(:,3)/Er)
ylabel('Er')
title('0.2 Er')

[qs,Es]=OneDBand(V2,lambda/2);
subplot(2,3,2)
plot(qs,Es(:,1)/Er)
hold on;
plot(qs,Es(:,2)/Er)
hold on;
plot(qs,Es(:,3)/Er)
ylabel('Er')
title('0.7 Er')

[qs,Es]=OneDBand(V3,lambda/2);
subplot(2,3,3)
plot(qs,Es(:,1)/Er)
hold on;
plot(qs,Es(:,2)/Er)
hold on;
plot(qs,Es(:,3)/Er)
ylabel('Er')
title('1.2 Er')

[qs,Es]=OneDBand(V4,lambda/2);
subplot(2,3,4)
plot(qs,Es(:,1)/Er)
hold on;
plot(qs,Es(:,2)/Er)
hold on;
plot(qs,Es(:,3)/Er)
ylabel('Er')
title('3 Er')

[qs,Es]=OneDBand(V5,lambda/2);
subplot(2,3,5)
plot(qs,Es(:,1)/Er)
hold on;
plot(qs,Es(:,2)/Er)
hold on;
plot(qs,Es(:,3)/Er)
ylabel('Er')
title('5 Er')

[qs,Es]=OneDBand(V6,lambda/2);
subplot(2,3,6)
plot(qs,Es(:,1)/Er)
hold on;
plot(qs,Es(:,2)/Er)
hold on;
plot(qs,Es(:,3)/Er)
ylabel('Er')
title('10 Er')