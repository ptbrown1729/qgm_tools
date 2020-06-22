function [Img] = transformObj2Img(Ximg, Yimg, P, Obj, Xobj, Yobj)
%[Img] = transformObj2Img(Ximg,Yimg,P,Obj,Xobj,Yobj)
%Ximg and Yimg are a grid of the points to evaluate Img at...so Img is the
%same size as these
%P = [Cx_i,Cy_i,Theta,Scale] are the parameters of the transformation
%Img is the image space representation, at coordinates Ximg, Yimg.
%if Xobj and Yobj are not specified, it is assumed that the center of the
%pictures is zero.
%%%General
%The transformation and inverse transformation are given by...
%[Xobj;Yobj] = [[cos(Theta) -sin(Theta)];[sin(Theta) cos(Theta)]]*([Ximg;Yimg]-[Cx_img;Cy_img])*1/Scale
%[Ximg;Yimg] = Scale*[[cos(Theta) sin(Theta)];[-sin(Theta) cos(Theta)]]*[Xobj;Yobj] + [Cx_img;Cy_img]

    if ~exist('Xobj','var')||~exist('Xobj','var')
        %if no coordinates provided for Obj, assume that it is centered
        %about the origin.
        x = (1:size(Obj, 2)) - (size(Obj, 2) + 1) / 2;
        y = (1:size(Obj, 1)) - (size(Obj, 1) + 1) / 2;
        [Xobj, Yobj] = meshgrid(x,y);
    end
    
    Cx_i = P(1);
    Cy_i = P(2);
    Theta = P(3);
    Scale = P(4);
    
%     Input coordinates transformed to image space...
%     Ximg = Scale*Xobj*cos(Theta) + Scale*Yimg*sin(Theta) + Cx_i;
%     Yimg = Scale*Yobj*cos(Theta) - Scale*Xobj*sin(Theta) + Cy_i;    
    C_o = 1/Scale * [[cos(Theta), -sin(Theta)]; [sin(Theta), cos(Theta)]] * [Cx_i; Cy_i];
    
    ObjInterp = @(Xo,Yo) interp2(Xobj, Yobj, Obj, Xo, Yo);
    ImgFn  = @(Xi,Yi) ObjInterp((Xi-Cx_i) * cos(Theta) / Scale - (Yi - Cy_i) * sin(Theta) / Scale,...
                                (Xi-Cx_i) * sin(Theta) / Scale + (Yi - Cy_i) * cos(Theta) / Scale);
    Img = ImgFn(Ximg, Yimg);
    Img(isnan(Img)) = 0;
    
end