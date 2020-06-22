w = 50:0.1:70;

dist = 11e-3;

s = w.*sqrt(1+(dist./(pi*(w*1e-6).^2/1064e-9)).^2);

plot(w,s)

%%
w = 60e-6;

dist = 100:1:400;

s = w.*sqrt(1+(dist*1e-3./(pi*w.^2/1064e-9)).^2)*1e3;

plot(dist,s)