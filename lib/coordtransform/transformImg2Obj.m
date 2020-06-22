function [Obj] = transformImg2Obj(Xobj, Yobj, P, Img, Ximg, Yimg)
%[Obj] = transformImg2Obj(Xobj,Yobj,P,Img,Ximg,Yimg)
%%%Args
%Xobj and Yobj are a grid of the points to evaluate Obj at...so Obj is the
%same size as these
%
%P = [Cx_i, Cy_i, Theta, Scale] are the parameters of the transformation
%
%Img is the image space representation, at coordinates Ximg, Yimg.
%if Ximg and Yimg are not specified, it is assumed that the center of the
%pictures is zero.
%
%%%General
%The transformation and inverse transformation are given by...
%[Xobj;Yobj] = [[cos(Theta) -sin(Theta)];[sin(Theta) cos(Theta)]]*([Ximg;Yimg]-[Cx_img;Cy_img])*1/Scale
%[Ximg;Yimg] = Scale*[[cos(Theta) sin(Theta)];[-sin(Theta) cos(Theta)]]*[Xobj;Yobj] + [Cx_img;Cy_img]

    if ~exist('Ximg','var')||~exist('Yimg','var')
        %if the coordinates of the image aren't provided, assume the image
        %is centered about the origin in image space.
        x = (1:size(Img, 2)) - (size(Img, 2) + 1) / 2;
        y = (1:size(Img, 1)) - (size(Img, 1) + 1) / 2;
        [Ximg, Yimg] = meshgrid(x, y);
    end
    
    Cx_i = P(1);
    Cy_i = P(2);
    Theta = P(3);
    Scale = P(4);
    
    ImgInterp = @(Xi, Yi) interp2(Ximg, Yimg, Img, Xi, Yi);
    ObjFn  = @(Xo, Yo) ImgInterp(Xo * cos(Theta) * Scale + Yo * sin(Theta) * Scale + Cx_i,...
                                -Xo * sin(Theta) * Scale + Yo * cos(Theta) * Scale + Cy_i);
 
%     Input coordinates transformed to object space...
%     Xobj = (Ximg-Cx_i)*cos(Theta)/Scale -(Yimg-Cy_i)*sin(Theta)/Scale;
%     Yobj = (Ximg-Cx_i)*sin(Theta)/Scale+(Yimg-Cy_i)*cos(Theta)/Scale;
    
    Obj = ObjFn(Xobj, Yobj);
    Obj(isnan(Obj)) = 0;

end