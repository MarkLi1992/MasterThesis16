classdef Solid8layered < handle
    
    properties
        elprop;
        lamIr
        ex; ey; ez;
        
        %Internal mesh properties
        nel
        intNno
        intNdofs
        intDof
        intEdof
        element2layer
        
        %Z coordinate for all elements, for example: [-1 -0.9 ... 1]
        %Used when computing stresses
        elementZcoords
        
        %Interpolator and Integrationrule for each element
        interp
        ir
    end
    
    methods
        
        function obj = Solid8layered(ex,ey,ez, elprop, nlamel)
            
            %Internal mesh
            nlay = elprop.nLam;
            obj.nel = nlay*nlamel;
            obj.intNno = obj.nel*4 + 4;
            obj.intNdofs = obj.intNno*3;
            obj.intDof = reshape(1:obj.intNdofs,3,obj.intNno)';
            obj.intEdof = zeros(24,obj.nel);
            temp = 1:4;
            for i=1:obj.nel
                elNodes = [temp + 4*(i-1), temp + 4*i]';
                elDofs = obj.intDof(elNodes,:)';
                elDofs = elDofs(:);
                obj.intEdof(:,i) = elDofs;
                
                obj.ex(:,i) = ex(elNodes);
                obj.ey(:,i) = ey(elNodes);
                obj.ez(:,i) = ez(elNodes);
            end
            
            obj.element2layer = zeros(1,nlay);
            itr = 1;
            for ie = 1:nlay
                for il= 1:nlamel
                    obj.element2layer(itr) = ie;
                    itr = itr+1;
                end
            end
            
            elementZcoords = ez(1:4:end)';
            obj.elementZcoords = (2*elementZcoords - (elementZcoords(end)+elementZcoords(1)))/(elementZcoords(end)-elementZcoords(1));
            
            %Always simplest ir and intetp
            obj.interp = InterpolatorX2Y2Z2;
            obj.ir = IntegrationRule;
            obj.ir.setupCubeRule(2,2,2);
            
            obj.elprop = elprop;
        end
        
        function [K,f] = computeKandf(obj,eq)
            
            K = zeros(obj.intNdofs);
            f = zeros(obj.intNdofs,1);
            
            for iel = 1:obj.nel
                clay = obj.element2layer(iel);
                cedof = obj.intEdof(:,iel);
                cD = obj.elprop.Dmatrices(:,:,clay);
                cex = obj.ex(:,iel); cey = obj.ey(:,iel); cez = obj.ez(:,iel);
                
                
%                 eq = [0,0,0]';
                [Kout, fout] = solid8LayeredElement(cex',cey',cez',cD,eq, obj.interp, obj.ir);
        
                K(cedof,cedof) = K(cedof,cedof) + Kout;
                f(cedof) = f(cedof) + fout;

            end
            
        end
        
        function [stresses, zCoord, xCoord, yCoord] = computeStressThroughThickness(obj,a, local_points)
            
            npoints = size(local_points,2);
            
            for ip = 1:npoints;
                stressItr = 1;
                for iel = 1:obj.nel
                    clay = obj.element2layer(iel);
                    cedof = obj.intEdof(:,iel);
                    cD = obj.elprop.Dmatrices(:,:,clay);
                    cex = obj.ex(:,iel); cey = obj.ey(:,iel); cez = obj.ez(:,iel);
                    
                    ca = a(cedof);
                    zetaCoord = linspace(-1,1,3);
                    for ii = 1:length(zetaCoord);
                        cp = [local_points(1,ip),local_points(2,ip),zetaCoord(ii)];
                        zCoord(stressItr) = obj.interp.eval_N(cp) * cez;
                        stresses(ip).stress(:,stressItr) = solid8layeredStress(cex,cey,cez,ca,cD,cp, obj.interp);
                        stressItr = stressItr + 1;
                    end

                end
                xCoord(ip) = obj.interp.eval_N(cp) * cex;
                yCoord(ip) = obj.interp.eval_N(cp) * cey;
            end

        end
        
        function [stresses, outcoords] = computeStressAt(obj,a,coords)
            warning('This function Assumes that one layer is used.')
            npoints = size(coords,2);
            
            for ip = 1:npoints;
                
                for iel = 1
                    clay = obj.element2layer(iel);
                    cedof = obj.intEdof(:,iel);
                    cD = obj.elprop.Dmatrices(:,:,clay);
                    cex = obj.ex(:,iel); cey = obj.ey(:,iel); cez = obj.ez(:,iel);
                    
                    ca = a(cedof);
                    
                    cp = coords(:,ip);
                    
                    outcoords(:,ip) = obj.interp.eval_N(cp) * [cex,cey,cez];
                    stresses(ip).stress = solid8layeredStress(cex,cey,cez, ca, cD, cp, obj.interp);
                end
            end
            
        end
        
        function [T, outcoords] = computeShearForceAtXY(obj,a,coords, ez)
            %Definition of coords: -->  coords = [x y];
            
            warning('This function Assumes that one layer is used.')
            interp1D = InterpolatorX2;
            ir1d = IntegrationRule;
            ngp_xi = 2;
            ir1d.setupLineRule( ngp_xi);
            
            T = 0;
            for gp = ir1d.gps
                
                lcoords = [coords gp.local_coords];
                [~, detJ] = interp1D.eval_dNdx(gp.local_coords, ez);
                [stresses, outcoords] = obj.computeStressAt(a,lcoords');
                
                dz = gp.weight * detJ;
                T = T + stresses(1).stress(5) * (dz);
            end
        end
        
    end
    
end

%function for calculating stress given a coordinate
function [stress] = solid8layeredStress(ex,ey,ez,a,D,coord, interp)

%Gauss points
xi=coord(1);
eta=coord(2);
zeta=coord(3);

%Shape functions
[dNdx, ~, JT] = interp.eval_dNdx(coord, ex', ey', ez');

[JT0] = interp.eval_ContraBaseVectors([0,0,0], ex', ey', ez');
% [JT0] = interp.eval_ContraBaseVectors(coord, ex', ey', ez');
T0 = transMat( JT0 );

%Bmatrix
B = solid8Bmatrix(dNdx);

% stress = D*B*a;
stress = (T0^-1)*D*B*a;
% stress = B*a;
end

%The standard solid-elemenet stiffness
function [ Ke,fe ] = solid8LayeredElement(ex,ey,ez,D,eq, interp, ir)

%Init vectors
Ke = zeros(24,24);
fe = zeros(24,1);

%Start loop
for gp= ir.gps
    
    %Shape functions
    Nxieta = interp.eval_N(gp.local_coords);
    [dNdx, detJ] = interp.eval_dNdx(gp.local_coords, ex, ey, ez);
    [N,B] = solid8NandBmatrix(Nxieta,dNdx);
    
    %derivative
    Ke = Ke + B'*D*B * (detJ*gp.weight);
    fe = fe + N'*eq * (detJ*gp.weight);
    
end

end


