classdef Solid8EasLayered < handle
    
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
        
        %Submatrices are stored in case of EAS is used, and alpha is needed
        submatrices;
        
        %EAS matrix
        M;
    end
    
    methods
        
        function obj = Solid8EasLayered(ex,ey,ez, elprop, nlamel, M)
            
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
            obj.ir.setupCubeRule(3,3,3);
            
            obj.elprop = elprop;
            obj.M = M;
        end
        
        function [K,f] = computeKandf(obj,eq)
            
            K = zeros(obj.intNdofs);
            f = zeros(obj.intNdofs,1);
            
            for iel = 1:obj.nel
                clay = obj.element2layer(iel);
                cedof = obj.intEdof(:,iel);
                cD = obj.elprop.Dmatrices(:,:,clay);
                cex = obj.ex(:,iel); cey = obj.ey(:,iel); cez = obj.ez(:,iel);
                cM = obj.M;
                
%                 eq = [0,0,0]';
                [Kout, fout, Ke, He, Le ] = solid8EasLayeredElement(cex',cey',cez',cD,cM,eq, obj.interp, obj.ir);
                
                K(cedof,cedof) = K(cedof,cedof) + Kout;
                f(cedof) = f(cedof) + fout;
                
                
                obj.submatrices(iel).Ke = Ke;
                obj.submatrices(iel).He = He;
                obj.submatrices(iel).Le = Le;
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
                    cM = obj.M;
                    
                    ca = a(cedof);
                    zetaCoord = linspace(-1,1,2);
                    for ii = 1:length(zetaCoord);
                        cp = [local_points(1,ip),local_points(2,ip),zetaCoord(ii)];
                        zCoord(stressItr) = obj.interp.eval_N(cp) * cez;
                        cAlpha = obj.computeAlpha(ca, iel);
                        stresses(ip).stress(:,stressItr) = solid8EasLayeredStress(cex,cey,cez, ca,cAlpha, cD, cM, cp, obj.interp);
                        stressItr = stressItr + 1;
                    end
                end
                xCoord(ip) = obj.interp.eval_N(cp) * cex;
                yCoord(ip) = obj.interp.eval_N(cp) * cey;
            end
        end
        
        function [alpha] = computeAlpha(obj, a, iel)
            %Clarification of the code below: -----> alpha = inv(He)*Le*a
            alpha = -inv(obj.submatrices(iel).He)*obj.submatrices(iel).Le*a;
        end
        
        function [stresses, outcoords] = computeStressAt(obj,a,coords)
            warning('This function Assumes that one layer is used.')
            npoints = size(coords,2);
            
            for ip = 1:npoints;
                
                %Copied from function "computeStressThroughThickness" so
                %there might me sum unnecessary code
                for iel = 1
                    clay = obj.element2layer(iel);
                    cedof = obj.intEdof(:,iel);
                    cD = obj.elprop.Dmatrices(:,:,clay);
                    cex = obj.ex(:,iel); cey = obj.ey(:,iel); cez = obj.ez(:,iel);
                    cM = obj.M;
                    
                    ca = a(cedof);
                    
                    cp = coords(:,ip);
                    
                    outcoords(:,ip) = obj.interp.eval_N(cp) * [cex,cey,cez];
                    cAlpha = obj.computeAlpha(ca, iel);
                    stresses(ip).stress = solid8EasLayeredStress(cex,cey,cez, ca,cAlpha, cD, cM, cp, obj.interp);
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
        
        function [M, outcoords] = computeMomentAtXY(obj,a,coords, ez)
            %Definition of coords: -->  coords = [x y];
            
            warning('This function Assumes that one layer is used.')
            interp1D = InterpolatorX2;
            ir1d = IntegrationRule;
            ngp_xi = 4;
            ir1d.setupLineRule( ngp_xi);
            
            M = 0;
            for gp = ir1d.gps
                
                lcoords = [coords gp.local_coords];
                [~, detJ] = interp1D.eval_dNdx(gp.local_coords, ez);
                [stresses, outcoords] = obj.computeStressAt(a,lcoords');
                
                dz = gp.weight * detJ;
                M = M + stresses(1).stress(1)* gp.local_coords * (dz);
            end
        end
        
        function [alpha] = computeEpsilonAlphaValue(obj, a, ez)

            interp1D = InterpolatorX2;
            ir1d = IntegrationRule;
            ngp_xi = 3;
            ir1d.setupLineRule( ngp_xi);
            
            Tforce = 0;
            for gp = ir1d.gps
                
                [~, detJ] = interp1D.eval_dNdx(gp.local_coords, ez);
                zcoord = interp1d.eval_N(gp.local_coords) * [obj.elprop.int_coordsG(1) obj.elprop.int_coordsG(2)]
                dz = gp.weight * detJ;
                Tforce = Tforce + obj.computeTauIntegral(zcoord) * (dz);
            end
            
            computeShearForceAtXY(obj,a, [0,0], ez)
            
            epsAlpha ;
            
        end
        
        function [Tau] = computeTauIntegral(obj, zcoord)
            
%             upToEl = getElementFromZCoord(zcoord);
            Tau = 0;
            zm = (obj.elprop.int_coordsG(1) + obj.elprop.int_coordsG(end))/2;
            %integrate elements from top to bottom
            for iel=1:obj.nel%upToEl-1
                clay = obj.element2layer(iel);
                cD = obj.elprop.Dmatrices(:,:,clay);
                ztop = obj.elprop.int_coordsG(iel+1); 
                zbot = obj.elprop.int_coordsG(iel);
                
                 %"Integral":
                if(zcoord>zbot && zcoord<ztop)
                    Tau = cD(1,1)*(0.5*(zcoord^2 - zbot^2) + zm*(zbot - zcoord)) + Tau;
                    break;
                else                 
                   Tau = cD(1,1)*(0.5*(ztop^2 - zbot^2) + zm*(zbot - ztop)) + Tau;
                end
  
            end

        end
    end
    
end

function [stress] = solid8EasLayeredStress(ex,ey,ez,a,alpha,D,Mi,coord, interp)

%Gauss points
xi=coord(1);
eta=coord(2);
zeta=coord(3);

%Shape functions
[dNdx, detJ] = interp.eval_dNdx(coord, ex', ey', ez');

%Bmatrix
B = solid8Bmatrix(dNdx);

%EAS part
[~, detJ0, JT0] = interp.eval_dNdx([0,0,0], ex', ey', ez');
Mtmp = Mi(xi,eta,zeta);
T0 = transMat( JT0' );
M = detJ0/detJ * T0 * Mtmp;

[JT0] = interp.eval_ContraBaseVectors([0,0,0], ex', ey', ez');
T = transMat( JT0 );

stress = D*B*a + D*M*alpha;
stress = (T^-1)*stress;
end

%The standard solid-elemenet with EAS
function [ Kout,fout, Ke, He, Le ] = solid8EasLayeredElement(ex,ey,ez,D,Mi, eq,interp, ir)

%Eas ndofs
nEnhDofs = size(Mi(0,0,0),2);

%Init vectors
Ke = zeros(24,24);
He = zeros(nEnhDofs,nEnhDofs);
Le = zeros(nEnhDofs,24);
fe = zeros(24,1);

[~, detJ0, JT0] = interp.eval_dNdx([0,0,0], ex, ey, ez);
T0 = transMat( JT0' );
%Start loop
for gp= ir.gps
    
    %Shape functions
    Nxieta = interp.eval_N(gp.local_coords);
    [dNdx, detJ, JT] = interp.eval_dNdx(gp.local_coords, ex, ey, ez);
    [N,B] = solid8NandBmatrix(Nxieta,dNdx);
    
    %Enanced part
    Mtemp = Mi(gp.local_coords(1),gp.local_coords(2),gp.local_coords(3));
    M = detJ0/detJ * T0*Mtemp;
    
    %Integrate
    Ke = Ke + B'*D*B * (detJ*gp.weight);
    He = He + M'*D*M * (detJ*gp.weight);
    Le = Le + M'*D*B * (detJ*gp.weight);
    
    fe = fe + N'*eq * (detJ*gp.weight);
    
end

Kout = Ke - Le'*inv(He)*Le;
fout = fe;

end



