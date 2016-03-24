classdef Mesh < handle
    %Mesh class to store all data associated with the mesh
    % 
    
    properties 
        lx; ly; lz;
        nelx; nely; nelz;

        edof
        nomesh %node mesh conectivity
        coord
        ex
        ey
        ez
        dof
        nel
        ndofs
        nno
        neldofs
        sideElements
        
        %Variables for stacked solid8 elements for composite modeling
        %These are not necessarily used
        nlam; nlamel;
    end
    
    methods
        function create_cube_mesh(obj,lx,ly,lz,nelx,nely,nelz)
            [obj.edof,... 
            obj.coord,...
            obj.ex,...
            obj.ey,...
            obj.ez,...
            obj.dof,...
            obj.nel,...
            obj.ndofs,...
            obj.nno,...
            obj.nomesh,...
            obj.sideElements] = ...
            cubeMesher(lx,ly,lz,nelx,nely,nelz);
            
            obj.neldofs = size(obj.edof,1); 
            obj.nelx = nelx; obj.nely = nely; obj.nelz = nelz;
            obj.lx = lx; obj.ly = ly; obj.lz = lz;

        end
        
        function create_cube_mesh_stacked_solid_elements(obj,lx,ly,lz,nelx,nely,nlam, nlamel)
            [obj.edof,... 
            obj.coord,...
            obj.ex,...
            obj.ey,...
            obj.ez,...
            obj.dof,...
            obj.nel,...
            obj.ndofs,...
            obj.nno,...
            obj.sideElements] = ...
            cubeMesherComposite(lx,ly,lz,nelx,nely,nlam,nlamel);
            
            obj.neldofs = size(obj.edof,1); 
            obj.nelx = nelx; obj.nely = nely; obj.nlam = nlam; obj.nlamel = nlamel;
            obj.lx = lx; obj.ly = ly; obj.lz = lz;

        end
        
    end
    
end

