classdef Solid8AnsEasSR < handle
    
    properties
        elprop;
        lamIr
        ex; ey; ez;
        
        %Internal mesh properties
        lx;ly;lz; 
        lx3;ly3;lz3;
        
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
        
        function obj = Solid8AnsEasSR(ex,ey,ez, elprop, M)
            obj.elprop = elprop;
            obj.M = M;
            
            ex3 = [ex(1:4);ex]; ey3 = [ey(1:4);ey];
            [obj.lx, obj.ly, obj.lz] = ex2lx(ex,ey,ez,obj.elprop.int_coordsG);
           
            for i=1:elprop.nLam
                 obj.lx3(:,i) = ex3;
                 obj.ly3(:,i) = ey3;
                 obj.lz3(:,i) = [[1 1 1 1]*elprop.int_coordsG(i), [1 1 1 1]*mean(elprop.int_coordsG([i, i+1])), [1 1 1 1]*elprop.int_coordsG(i+1)]';
            end
            obj.ex = ex; obj.ey = ey; obj.ez = ez;
            
            %Always simplest ir and intetp
            obj.interp = InterpolatorX2Y2Z2;
            obj.ir = IntegrationRule;
            obj.ir.setupCubeRule(3,3,3);
            
            obj.lamIr = LayeredIntegrationRule(elprop.nLam, 3,3,3);
            
            nEnhDofs = size(M(0,0,0),2);
            obj.submatrices.Ke = zeros(24);
            obj.submatrices.He = zeros(nEnhDofs);
            obj.submatrices.Le = zeros(nEnhDofs,24);
        end
       
        function [K,f] = computeKandf(obj, eq)
            
            nEnhDofs = size(obj.M(0,0,0),2);
            interpNormal = InterpolatorX2Y2;
            
            %Init vectors
            Ke = zeros(24,24);
            He = zeros(nEnhDofs,nEnhDofs);
            Le = zeros(nEnhDofs,24);
            fe = zeros(24,1);
            
            [~, detJ0, JT0] = obj.interp.eval_dNdx([0,0,0], obj.ex', obj.ey', obj.ez');
            T0 = transMat( JT0' );
                                
            %Interpolation points for ANS
            [B13a,B23b,B13c,B23d] = BAnsShear(obj.interp,  obj.ex', obj.ey', obj.ez', inv(T0));
            [B33e,B33f,B33g,B33h] = BAnsNormal(obj.interp, obj.ex', obj.ey' ,obj.ez', inv(T0));
            
            for ilay = 1:obj.elprop.nLam
                D = obj.elprop.Dmatrices(:,:,ilay);
                for gp = obj.lamIr.irs(ilay).gps
                    %Get coords
                    lcoords = gp.local_coords;
                    xi = lcoords(1); eta = lcoords(2); zeta = lcoords(3);
                    
                    %Transform gp-coord to element-coordinate
                    ecoords = lcoords;
                    ecoords(3) = obj.lamIr.getElementGaussCoordinate(lcoords(3), obj.elprop.int_coordsL(ilay:ilay+1));
                    
                    %Shape functions
                    Nxieta = obj.interp.eval_N(ecoords);
                    [dNdx, ~] = obj.interp.eval_dNdx(ecoords, obj.ex', obj.ey', obj.ez');
                    [N,~] = solid8NandBmatrix(Nxieta,dNdx);
                    
                    %Interpolation matrix for ANS   
                    Btmp = inv(T0)*solid8Bmatrix(dNdx);
                    
                    B13tilde = 0.5*(1-eta)*B13a + 0.5*(1+eta)*B13c; 
                    Btmp(5,:) = B13tilde;

                    B23tilde = 0.5*(1-xi)*B23d + 0.5*(1+xi)*B23b; 
                    Btmp(6,:) = B23tilde;

                    N33 = interpNormal.eval_N(gp.local_coords([1 2]));
                    B33tilde =  N33(1)*B33e + N33(2)*B33f + N33(3)*B33h + N33(4)*B33g;
                    Btmp(3,:) = B33tilde;
                    B = T0*Btmp;
                    
                    %Determinant of j for integration
                    [~, detJ] = obj.interp.eval_dNdx(ecoords, obj.lx(:,ilay)', obj.ly(:,ilay)', obj.lz(:,ilay)');

                    %Enanced part
                    Mtemp = obj.M(ecoords(1),ecoords(2),ecoords(3));
                    Mi = detJ0/detJ * T0*Mtemp;
                    
                    %Integrate
                    Ke = Ke + B'*D*B * (detJ*gp.weight);
                    He = He + Mi'*D*Mi * (detJ*gp.weight);
                    Le = Le + Mi'*D*B * (detJ*gp.weight);
                    
                    fe = fe + N'*eq * (detJ*gp.weight);
                end
            end
            K = Ke - Le'*inv(He)*Le;
            f = fe;
            
            obj.submatrices.Ke = Ke;
            obj.submatrices.Le = Le;
            obj.submatrices.He = He;
            
        end
           
        function [stresses, zcoord, plotStress, plotZ] = computeStressThroughLayers(obj,a, local_point)
            %input:
            % local_points = [x1 x2 ...
            %                 y1 y2 ...];
            alpha = obj.computeAlpha(a);
            stresses = zeros(6,2,obj.elprop.nLam);
            for ilay=1:obj.elprop.nLam
                
                %Loop over thickness points;
                z1 = obj.elprop.int_coordsL(ilay);
                z2 = obj.elprop.int_coordsL(ilay+1);

                coord1 =  [ local_point z1];
                coord2 =  [ local_point z2];
                
                [stresses_bot] = solid8anseasStress(obj.ex,obj.ey,obj. ez,a,alpha,...
                                                       obj.elprop.Dmatrices(:,:,ilay), obj.M, coord1, obj.interp);
                [stresses_top] = solid8anseasStress(obj.ex, obj.ey, obj.ez, a,alpha,...
                                                       obj.elprop.Dmatrices(:,:,ilay), obj.M, coord2, obj.interp);
                stresses(:,:,ilay) = [stresses_bot, stresses_top];
                int = InterpolatorX2Y2Z2;
                
                plotStress( :, (1:2) + (2*(ilay-1)) ) = [stresses_bot, stresses_top];
                plotZ( (1:2) + (2*(ilay-1)) ) = obj.elprop.int_coordsG( (1:2) + (ilay-1));
%                 error('output u as')
            end
            zcoord = obj.elprop.int_coordsG;
        end
        
        function [stresses, outcoords] = computeStressThroughThickness(obj,a, local_points)
            %input:
            % local_points = [x1 x2 ...
            %                 y1 y2 ...];
            
            npoints = size(local_points,2);
            
            %Loop over the points 
            for ip=1:npoints
                iitr = 1;
                %Loop over thickness points;
                zz = linspace(-1,1,100);
                for iz = 1:length(zz)
                    currentCoord =  [ local_points(:,ip); zz(iz)];
                    [stresses(ip).stress(:,iitr), ~, outcoords(:,iitr)] = obj.computeStressAt(a,currentCoord);
                    iitr = iitr + 1;
                end
                
            end
            
            
        end
        
        function [alpha] = computeAlpha(obj, a)
            %Clarification of the code below: -----> alpha = inv(He)*Le*a
            alpha = -inv(obj.submatrices.He)*obj.submatrices.Le*a;
        end
        
        function [stresses, strains, outcoords] = computeStressAt(obj,a,coords)
            
            %Get current layer at coordinates in "coords"
            clay = obj.getLayerFromZCoord(coords(3),'local');

            %D-matrix for this layer
            cD = obj.elprop.Dmatrices(:,:,clay);
            
            %Global coordinates
            outcoords = obj.interp.eval_N(coords) * [obj.ex,obj.ey,obj.ez];
            cAlpha = obj.computeAlpha(a);
%             [stresses, strains] = solid8EasLayeredStress(obj.ex,obj.ey,obj.ez, a, cAlpha, cD, obj.M, coords, obj.interp);
            [stresses] = solid8anseasStress    (obj.ex,obj.ey,obj.ez, a, cAlpha, cD, obj.M, coords, obj.interp);
            strains = stresses*0;
        end 
        
        function [] = postProcess2(obj, as)
            %Calculate witdth and length of element
            elx = abs((obj.ex(1) - obj.ex(2)));
            ely = abs((obj.ey(1) - obj.ey(3)));
            
            tauItr = 1;
            Taubc = [0,0]';
            Taubc = [0]';
            for ilay = 1:obj.elprop.nLam
                DLT = obj.elprop.Dmatrices([1,2,4],[1 2 4], ilay)
                
                zminus = obj.elprop.int_coordsG(ilay);
                zplus  = obj.elprop.int_coordsG(ilay+1);
                layer_zminus = obj.elprop.int_coordsL(ilay);
                layer_zplus  = obj.elprop.int_coordsL(ilay+1);
                
                keyboard
                local_zpoints = linspace(layer_zminus,layer_zplus,100); %could map the points instread but this is easier
                zpoints = linspace(zminus,zplus,100);

                for iz = 2:length(zpoints);
                                    
                    %Integrate the equilbrium equation
                    interp1d = InterpolatorX2; ir1d = IntegrationRule; ir1d.setupLineRule(3); integral_temp = 0;
                    for gp=ir1d.gps
                        [~,detJ] = interp1d.eval_dNdx(gp.local_coords, [zminus, zpoints(iz)]);
                        stressEvalPoint = obj.lamIr.getElementGaussCoordinate(gp.local_coords, [layer_zminus, local_zpoints(iz)]');
%                         stressGrads = evalStressGradients_v2(obj.interp, [stressEvalPoint], as, DLT, elx,ely);
                        stressGrads = evalStressGradientsFromStrain(obj.interp, [0,0,stressEvalPoint], as, DLT, obj.ex, obj.ey, obj.ez);
                        dz = detJ * gp.weight;
                        integral_temp = integral_temp + (stressGrads.sigxx_x + stressGrads.sigxy_x) * dz;
                    end
                    
                    stressEvalPoint = obj.lamIr.getElementGaussCoordinate(0, [layer_zminus, local_zpoints(iz)]');
                    stressGrads = evalStressGradientsFromStrain(obj.interp, [0,0,stressEvalPoint], as, DLT, obj.ex, obj.ey, obj.ez);                 
                    stressgrad_save(tauItr) = stressGrads.sigxx_x;
                    
                    sig_xz(1, tauItr) = integral_temp + Taubc;
                    
                    zplot(1, tauItr) = zpoints(iz);
                    tauItr = tauItr+1;
                end
                
                Taubc = sig_xz(:,end);
            end
            figure
            plot(sig_xz, zplot);
            
            figure;
            plot(stressgrad_save, zplot)
            keyboard;
        end
         
        function [tauxz, tauyz] = ShearStressesPostProcess(obj, as, opt)
            
            if(strcmp(opt,'shear'))
                dim = 3;
                stressInterpolator = InterpolatorX2Y2Z2;
            elseif(strcmp(opt,'normal'))
                stressInterpolator = InterpolatorX2Y2Z3;
                dim = 2;
            end
%             keyboard
            %Dim is the number of dofs per node. ie, if you want to post
            %process shear stresses, then dim = 3, and postprocess normal
            %stress dim = 2;
            
            %variable as is stresses in each node, for each layer
            tauItr = 1;
            if(strcmp(opt,'shear')); Taubc = [0,0]'; integral_temp = [0,0]'; elseif(strcmp(opt,'normal')); Taubc = 0; integral_temp = 0; end
            for ilay = 1:obj.elprop.nLam
                %Get properties for this layer
                current_as = as(:,ilay);
                
                %Get coordinates for this layer
                zminus = obj.elprop.int_coordsG(ilay);
                zplus  = obj.elprop.int_coordsG(ilay+1);
                layer_zminus = obj.elprop.int_coordsL(ilay);
                layer_zplus  = obj.elprop.int_coordsL(ilay+1);
                
                %Points where to evaluate the stress (in current layer)
                local_zpoints = linspace(layer_zminus,layer_zplus,3); %could map the points instread but this is easier
                zpoints = linspace(zminus,zplus,3);

                %Loop through the points
                for iz = 1:length(zpoints);
                                    
                    %Integrate the equilbrium equation
                    interp1d = InterpolatorX2; ir1d = IntegrationRule; ir1d.setupLineRule(3);
                    for gp=ir1d.gps
                        [~,detJ] = interp1d.eval_dNdx(gp.local_coords, [zminus, zpoints(iz)]);
                        stressEvalPoint = obj.lamIr.getElementGaussCoordinate(gp.local_coords, [layer_zminus, local_zpoints(iz)]');
                        stressEvalPoint = obj.lamIr.getLayerGaussCoordinate(stressEvalPoint, [layer_zminus, layer_zplus]');
                        
                        stressGrads = evalGradients(obj.lx(:,ilay), obj.ly(:,ilay), obj.lz(:,ilay), obj.interp, current_as, [0,0,stressEvalPoint] , dim);
                        dz = detJ * gp.weight;
                        
                        if(dim == 3) %If integration shear stress
                        integral_temp(1) = integral_temp(1) + -(stressGrads(1) + stressGrads(6)) * dz;
                        integral_temp(2) = integral_temp(2) + -(stressGrads(5) + stressGrads(3)) * dz;
                        elseif(dim == 2) %If integrating normal stress
                        integral_temp = integral_temp + (stressGrads(1) + stressGrads(3)) * dz;    
                        end
                        
%                         if(iz == 1) %Will produce value NaN.
%                             integral_temp(:) = 0;
%                         end
                    end
                    
                    %Save values
                    integral_save(:, tauItr) = integral_temp + Taubc;
                    zplot(1, tauItr) = zpoints(iz);
                    tauItr = tauItr+1;
                    integral_temp(:) = 0;
                end
                
                %Store boundary-stress for the next layer
                Taubc = integral_save(:,end);
            end
            integral_save(:,4:3:end) = [];
            integral_save = integral_save';
            tauxz = integral_save(:,1);
            tauyz = integral_save(:,2);
%             figure
%             plot(tauxz, zplot);
%             keyboard;
        end
        
        function [tauxz, tauyz] = ShearStressesPostProcess2(obj, as)
           
            dim = 6;
            
            %variable as is stresses in each node, for each layer
            tauItr = 1;
            Taubc = [0,0]'; 
            integral_temp = [0,0]'; 
            for ilay = 1:obj.elprop.nLam
                %Get properties for this layer
                current_as = as(:,ilay);
                
                %Get coordinates for this layer
                zminus = obj.elprop.int_coordsG(ilay);
                zplus  = obj.elprop.int_coordsG(ilay+1);
                layer_zminus = obj.elprop.int_coordsL(ilay);
                layer_zplus  = obj.elprop.int_coordsL(ilay+1);
                
                %Points where to evaluate the stress (in current layer)
                local_zpoints = linspace(layer_zminus,layer_zplus,3); %could map the points instread but this is easier
                zpoints = linspace(zminus,zplus,3);

                %Loop through the points
                for iz = 1:length(zpoints);
                                    
                    %Integrate the equilbrium equation
                    interp1d = InterpolatorX2; ir1d = IntegrationRule; ir1d.setupLineRule(3);
                    for gp=ir1d.gps
                        [~,detJ] = interp1d.eval_dNdx(gp.local_coords, [zminus, zpoints(iz)]);
                        stressEvalPoint = obj.lamIr.getElementGaussCoordinate(gp.local_coords, [layer_zminus, local_zpoints(iz)]');
                        stressEvalPoint = obj.lamIr.getLayerGaussCoordinate(stressEvalPoint, [layer_zminus, layer_zplus]');
                        
                        stressGrads = evalGradients2(obj.lx(:,ilay), obj.ly(:,ilay), obj.lz(:,ilay), obj.interp, current_as, [0,0,stressEvalPoint] , dim);
                        dz = detJ * gp.weight;
                        
                        integral_temp(1) = integral_temp(1) + -(stressGrads(1) + stressGrads(10)) * dz;
                        integral_temp(2) = integral_temp(2) + -(stressGrads(8) + stressGrads(4)) * dz;
                    end
                    
                    
                    %--------------
                    tmp = evalGradients2(obj.lx(:,ilay), obj.ly(:,ilay), obj.lz(:,ilay), obj.interp, current_as, [0 0 local_zpoints(iz)], dim);
                    syy_save(tauItr) = tmp(8);
                    %-
                    int = InterpolatorX2Y2Z2;
                    [JT0] = int.eval_ContraBaseVectors([0,0,local_zpoints(iz)], obj.lx(:,ilay)', obj.ly(:,ilay)', obj.lz(:,ilay)');
                    T = transMat( JT0 );
                    for i=1:8
                        ind = (1:6) + (6*(i-1));
                       newAs(ind,1) = T*current_as(ind,1);
                    end
                    current_as2 = newAs;
                    sigtmp = solid8NMatrix(int.eval_N([0,0,local_zpoints(iz)]),dim) * current_as2;
                    sigtmp_save(tauItr) = sigtmp(2);
                    %------------
                    
                    
                    %Save values
                    integral_save(:, tauItr) = integral_temp + Taubc;
                    zplot(1, tauItr) = zpoints(iz);
                    tauItr = tauItr+1;
                    integral_temp(:) = 0;
                end
                
                %Store boundary-stress for the next layer
                Taubc = integral_save(:,end);
            end
            integral_save(:,4:3:end) = [];
            integral_save = integral_save';
            tauxz = integral_save(:,1);
            tauyz = integral_save(:,2);
            
            
            %------
%             figure
%             sigtmp_save(4:3:end) = [];
%             zplot(4:3:end) = [];
%             plot(tauxz, zplot);
%             hold on; 
%             plot(tauyz, zplot);
%             
%             figure('name','Syy'); 
%             plot(sigtmp_save,zplot);
%             
%             
%             figure('name','Syy,y')
%             syy_save(4:3:end) = [];
%             plot(syy_save,zplot);
%             keyboard;
            %----
        end
        
        function [sig_zz] = NormalStressPostProcess(obj, as)
            
            TauInterp = InterpolatorX2Y2Z3;
            %variable as is stresses in each node, for each layer
            Itr = 1;
            bc = 0; 
            integral_temp = 0;
            for ilay = 1:obj.elprop.nLam
                %Get properties for this layer
                current_as = as(:,ilay);
                
                %Get coordinates for this layer
                zminus = obj.elprop.int_coordsG(ilay);
                zplus  = obj.elprop.int_coordsG(ilay+1);
                layer_zminus = obj.elprop.int_coordsL(ilay);
                layer_zplus  = obj.elprop.int_coordsL(ilay+1);
                
                %Points where to evaluate the stress (in current layer)
                local_zpoints = linspace(layer_zminus,layer_zplus,4); 
                zpoints = linspace(zminus,zplus,4);

                %Loop through the points
                for iz = 1:length(zpoints);
                                    
                    %Integrate the equilbrium equation
                    interp1d = InterpolatorX2; ir1d = IntegrationRule; ir1d.setupLineRule(5);
                    for gp=ir1d.gps
                        [~,detJ] = interp1d.eval_dNdx(gp.local_coords, [zminus, zpoints(iz)]);
                        stressEvalPoint = obj.lamIr.getElementGaussCoordinate(gp.local_coords, [layer_zminus, local_zpoints(iz)]');
                        stressEvalPoint = obj.lamIr.getLayerGaussCoordinate(stressEvalPoint, [layer_zminus, layer_zplus]');
                        
                       
                        stressGrads = evalGradients( obj.lx3(:,ilay), obj.ly3(:,ilay), obj.lz3(:,ilay), TauInterp, current_as, [0,0,stressEvalPoint] , 2);
                       
                        dz = detJ * gp.weight;                        
                        integral_temp = integral_temp + -(stressGrads(1) + stressGrads(3)) * dz;    

                    end
                    
                    %Save values
                    integral_save(:, Itr) = integral_temp + bc;
                    zplot(1, Itr) = zpoints(iz);
                    Itr = Itr+1;
                    integral_temp(:) = 0;
                end
                
                %Store boundary-stress for the next layer
                bc = integral_save(:,end);
            end
            integral_save(:,5:4:end) = [];
            sig_zz = integral_save';
%             keyboard;
%             figure
%             plot(integral_save(1,:), zplot);
%             keyboard;
        end
        
        function [el] = getLayerFromZCoord(obj,zc,opt)
            
            if(strcmp(opt,'global'))
              zcoords = obj.elprop.int_coordsG;
            elseif(strcmp(opt,'local'))
              zcoords = obj.elprop.int_coordsL;
            end
            
            for iz = 1:length(zcoords)-1
               
                if(zc >= zcoords(iz) && zc <= zcoords(iz+1))
                    el = iz;
                    return;
                end
                
            end
            
        end
        
    end
    
end

function stressgrads = evalStressGradientsFromStrain(interp, coord, as, DLT, lx, ly, lz)
    %Stress = Strains hehe
    
    %sigma xx
    dN = interp.eval_dNdx(coord, lx',ly',lz');
    dNdx = solid8NMatrix(dN(1,:),3);
    dNdy = solid8NMatrix(dN(2,:),3);
    
    sig_x = DLT*dNdx*as;
    sig_y = DLT*dNdy*as;
    
    stressgrads.sigxx_x = sig_x(1);
    stressgrads.sigyy_y = sig_y(2);
    stressgrads.sigxy_y = sig_y(3);
    stressgrads.sigxy_x = sig_x(3);
end

function stressgrads = evalGradients(lx, ly, lz, interp, as, coord, dim)
    %Derivatives of N-matrix
    dN = interp.eval_dNdx(coord, lx',ly',lz');
    dNdx = solid8NMatrix(dN(1,:),dim);
    dNdy = solid8NMatrix(dN(2,:),dim);
    
    sig_x = dNdx*as;
    sig_y = dNdy*as;
    
    stressgrads = [sig_x; sig_y];
%     stressgrads.sigxx_x = sig_x(1);
%     stressgrads.sigyy_y = sig_y(2);
%     stressgrads.sigxy_y = sig_y(3);
%     stressgrads.sigxy_x = sig_x(3);
end

function stressgrads = evalGradients2(lx, ly, lz, interp, as, coord, dim)
    
    [JT0] = interp.eval_ContraBaseVectors(coord, lx', ly', lz');
    T = transMat( JT0 );
% 
%     for i=1:8
%         ind = (1:6) + (6*(i-1));
%        newAs(ind,1) = T*as(ind,1);
%     end
%     as = newAs;

    %Derivatives of N-matrix
    [dNxy,~,Jscale] = interp.eval_dNdx(coord, lx',ly',lz');
%     dN = (JT0)*dN;
    dN = inv(JT0)*dNxy;%inv(Jscale)*interp.eval_dNdxi(coord);
    dNdx = solid8NMatrix(dN(1,:),dim);
    dNdy = solid8NMatrix(dN(2,:),dim);
    
    sig_x = dNdx*as;
    sig_y = dNdy*as;
    
    stressgrads = [sig_x; sig_y];
%     stressgrads.sigxx_x = sig_x(1);
%     stressgrads.sigyy_y = sig_y(2);
%     stressgrads.sigxy_y = sig_y(3);
%     stressgrads.sigxy_x = sig_x(3);
end

function stressgrads = evalStressGradients_v2(interp, zcoord, as, DLT, elx,ely)
    %Stress = Strains hehe
    
    %sigma xx
    sigxx1 = DLT*solid8NMatrix(interp.eval_N([-1 0 zcoord]),3) * as;
    sigxx2 = DLT*solid8NMatrix(interp.eval_N([ 1 0 zcoord]),3) * as;
    
    sigxx_x = (sigxx1(1) - sigxx2(1))/elx;
    sigxy_x = (sigxx1(3) - sigxx2(3))/elx;
    
    %sigma yy
    sigyy1 = DLT*solid8NMatrix(interp.eval_N([0 -1 zcoord]),3) * as;
    sigyy2 = DLT*solid8NMatrix(interp.eval_N([0  1 zcoord]),3) * as;
    
    sigyy_y = (sigyy1(2) - sigyy2(2))/ely;
    sigxy_y = (sigyy1(3) - sigyy2(3))/elx;
    
    stressgrads.sigxx_x = sigxx_x;
    stressgrads.sigyy_y = sigyy_y;
    stressgrads.sigxy_y = sigxy_y;
    stressgrads.sigxy_x = sigxy_x;
    
end

