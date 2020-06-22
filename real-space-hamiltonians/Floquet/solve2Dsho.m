[v,E]=floquetStep(@Hsho2D,1/6e3);

NN=13;

for ii=1:NN

subplot(ceil(NN/6),6,ii)
imagesc(real(reshape(v(:,end-ii),12,11)),[0,1])
axis equal;
end