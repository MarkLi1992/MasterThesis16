classdef Solid8AnsEasLayered < handle
    
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
        
        function obj = Solid8AnsEasLayered(ex,ey,ez, elprop, nlamel, M)
            
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
%                 [Kout, fout, Ke, He, Le ] = solid8AnsEasLayeredElement(cex',cey',cez',cD,cM,eq, obj.interp, obj.ir);
                [Kout, fout, Ke, fe, He, Le  ] = solid8anseas(cex,cey,cez, cD, cM, eq, [3,3,3]);
                
                K(cedof,cedof) = K(cedof,cedof) + Kout;
                f(cedof) = f(cedof) + fout;
                
                
                obj.submatrices(iel).Ke = Ke;
                obj.submatrices(iel).He = He;
                obj.submatrices(iel).Le = Le;
            end
            
        end
        
        function [stresses, zCoord, xCoord] = computeStressThroughThickness(obj,a, local_points)
            
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
                    zetaCoord = linspace(-1,1,10);
                    for ii = 1:length(zetaCoord);
                        cp = [local_points(1,ip),local_points(2,ip),zetaCoord(ii)];
                        zCoord(stressItr) = obj.interp.eval_N(cp) * cez;
                        cAlpha = obj.computeAlpha(ca, iel);
%                        stresses(ip).stress(:,stressItr) = solid8EasLayeredStress(cex,cey,cez, ca,cAlpha, cD, cM, cp, obj.interp);
                        stresses(ip).stress(:,stressItr) = solid8anseasStress(cex,cey,cez, ca,cAlpha, cD, cM, cp, obj.interp);
                        stressItr = stressItr + 1;
                    end
                end
                xCoord(ip) = obj.interp.eval_N(cp) * cex;
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

    end
    
end



