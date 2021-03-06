function Mesh = readMeshNC(filename,varargin)
% READMESHNC - Read a hydrodynamic model mesh file in NetCDF format.
% uvnode and nodexy are read in in whatever coordinate system they are in,
% but there are options to read ihn specifically named additional fields.
% 
% Inputs:   filename - filename of the NetCDF file (should be a ".nc" file)
%           varargin - Name-value pairs, options are
%               plotMesh - logical, plot the mesh?
%               readOS - logical, read OS grid refs (metres, not degrees)
%               readDeg - logical, read longitude and latitude into
%                       separate variable
%
% Outputs:  Mesh - a Matlab structure containing hydrodynamic mesh
%               variables
%
    plotMesh=0;
    readOS=0;
    readDeg=0;
    type=1;
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'plotMesh'
                plotMesh = varargin{i+1};
            case 'readOS'
                readOS = varargin{i+1};
            case 'readDeg'
                readDeg = varargin{i+1};
            case 'type'
                type = varargin{i+1};
        end
    end
    
    if (type==1)
        Mesh.uvnode=ncread(filename,'uvnode');
        Mesh.nodexy=ncread(filename,'nodexy');
        if (readOS==1)
            Mesh.uvnode_os=ncread(filename,'uvnode_os');
            Mesh.nodexy_os=ncread(filename,'nodexy_os');
        end
        if (readDeg==1)
            Mesh.uvnode_deg=ncread(filename,'uvnode_deg');
            Mesh.nodexy_deg=ncread(filename,'nodexy_deg');
        end
    
        
        Mesh.trinodes=ncread(filename,'trinodes');

        Mesh.uvdepth=ncread(filename,'depthUvnode');
        Mesh.depth=ncread(filename,'depthNodexy');
        
        Mesh.boundaryNodesAll=ncread(filename,'boundaryNodesAll');
        Mesh.boundaryNodesOpen=ncread(filename,'boundaryNodesOpen');
        Mesh.obc=ncread(filename,'boundaryNodesOpen');
        Mesh.siglay=ncread(filename,'siglay');
        
    elseif (type==2)
        lonc = ncread(filename,'lonc');
        latc = ncread(filename,'latc');
        Mesh.uvnode=[lonc,latc];
        
        lon = ncread(filename,'lon');
        lat = ncread(filename,'lat');
        Mesh.nodexy=[lon,lat];
        
        Mesh.trinodes=ncread(filename,'nv');
        Mesh.uvdepth=ncread(filename,'h_center');
        Mesh.depth=ncread(filename,'h');
        
        sl = ncread(filename,'siglay');
        Mesh.siglay=sl(1,:);
    end
    
    Mesh.nbe=ncread(filename,'nbe');

    if plotMesh
        plotMeshPDens(Mesh);
    end
end