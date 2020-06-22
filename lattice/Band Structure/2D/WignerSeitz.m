function vals = WignerSeitz(xpts,ypts,v1,v2)
%WignerSeitz, produces the characteristic function of the wigner seitz cell
%for a 2D lattice with basis vectors v1,v2.

xpts_padded=xpts;
ypts_padded=ypts;

if length(xpts)~=length(ypts)
    if (length(xpts)~=1)&&(length(ypts)~=1)
    	error('WignerSeitz:args','Problem with xpts and ypts arguments.');
    elseif length(xpts)==1
        xpts_padded=kron(xpts,ones(1,length(ypts)));
        ypts_padded=ypts;
    else
        xpts_padded=xpts;
        ypts_padded=kron(ypts,ones(1,length(xpts)));
    end
end


%could try to vectorize by making cell array and using cellfun...
vals=[];
for ii=1:length(xpts_padded)
    pt=[xpts_padded(ii),ypts_padded(ii)];
    testpts=[norm(pt-0),norm(pt-v1),norm(pt-v2),norm(pt+v1),norm(pt+v2),norm(pt-(v1+v2)),norm(pt-v1+v2),norm(pt-v2+v1),norm(pt+(v1+v2))];
    val=(norm(pt-0)==min(testpts));
    vals=cat(2,vals,val);
end


end

