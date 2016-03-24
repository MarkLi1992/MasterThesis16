function [ N ] = solid8NMatrix( Nxieta,dim )

if( ~(exist('dim','var')) ); dim = 3; end;

nshape = size(Nxieta,2);

for in=1:nshape
    N(:,(1:dim) + dim*(in-1)) = Nxieta(in)*eye(dim);
end

end