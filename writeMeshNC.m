function writeMeshNC(Mesh,filename,varargin)
    
    coordRef='WGS84';
    writeOS=0;
    obcFile=[];
    switchSigvec=0;
    %obcFile='C:\Users\SA01TA\Documents\MH_KTP\work\fvcomMeshData\WestCOMS\obc_nodes_minch.mat';
    
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'coordRef'
                coordRef = varargin{i+1};
            case 'writeOS'
                writeOS = varargin{i+1};
            case 'obcFile'
                obcFile = varargin{i+1};
            case 'switchSigvec' % reverses the sign of sigvec/siglay i.e. make negative if presently positive
                switchSigvec = varargin{i+1};
        end
    end
    
    %writeOS
    
    % Need to have SEPA's OS toolbox available on your path to do the
    % coordinate tansformation stuff
    %addpath(genpath('C:\Users\sa01ta\OneDrive - SAMS\Documents\code\matlab\os_toolbox-master'));
    
    % This bit ensures that:
    % - nodexy, uvnode are both in degrees in the file produced
    % - nodexy_os, uvnode_os are both in metres in the file produced
    if (strcmpi(coordRef,'OSGB1936'))
        [Mesh.nodexy_deg(:,1),Mesh.nodexy_deg(:,2)]=OS.convertAndTransform(Mesh.nodexy(:,1),Mesh.nodexy(:,2));
        [Mesh.uvnode_deg(:,1),Mesh.uvnode_deg(:,2)]=OS.convertAndTransform(Mesh.uvnode(:,1),Mesh.uvnode(:,2));
        Mesh.nodexy_os=Mesh.nodexy;
        Mesh.uvnode_os=Mesh.uvnode;
        Mesh.nodexy=Mesh.nodexy_deg;
        Mesh.uvnode=Mesh.uvnode_deg;
    else
        %disp('hello')
        Mesh.nodexy_deg=Mesh.nodexy;
        Mesh.uvnode_deg=Mesh.uvnode;
        if (writeOS==1)
            %disp('hello')
            [Mesh.nodexy_os(:,1),Mesh.nodexy_os(:,2)]=OS.convertAndTransform(Mesh.nodexy(:,1),Mesh.nodexy(:,2));
            [Mesh.uvnode_os(:,1),Mesh.uvnode_os(:,2)]=OS.convertAndTransform(Mesh.uvnode(:,1),Mesh.uvnode(:,2));
        end
    end

    delete(filename)
    % uvnode
    nccreate(filename,'uvnode','Dimensions',{'elems' length(Mesh.uvnode) 'two' 2},'Datatype','single');
    ncwrite(filename,'uvnode',Mesh.uvnode)
    % nodexy
    nccreate(filename,'nodexy','Dimensions',{'node' length(Mesh.nodexy) 'two' 2},'Datatype','single');
    ncwrite(filename,'nodexy',Mesh.nodexy)
    % uvnode_deg
    nccreate(filename,'uvnode_deg','Dimensions',{'elems' length(Mesh.uvnode_deg) 'two' 2},'Datatype','single');
    ncwrite(filename,'uvnode_deg',Mesh.uvnode_deg)
    % nodexy_deg
    nccreate(filename,'nodexy_deg','Dimensions',{'node' length(Mesh.nodexy_deg) 'two' 2},'Datatype','single');
    ncwrite(filename,'nodexy_deg',Mesh.nodexy_deg)
    if (writeOS==1)
        % % uvnode_os
        disp('hello')
        [E,N]=OS.convertAndTransform(Mesh.uvnode(:,1),Mesh.uvnode(:,2));
        nccreate(filename,'uvnode_os','Dimensions',{'elems' length(Mesh.uvnode) 'two' 2},'Datatype','single');
        ncwrite(filename,'uvnode_os',[E,N])
        % % nodexy_os
        [E,N]=OS.convertAndTransform(Mesh.nodexy(:,1),Mesh.nodexy(:,2));
        nccreate(filename,'nodexy_os','Dimensions',{'node' length(Mesh.nodexy) 'two' 2},'Datatype','single');
        ncwrite(filename,'nodexy_os',[E,N])
    end
    % trinodes
    nccreate(filename,'trinodes','Dimensions',{'nele' length(Mesh.uvnode) 'three' 3},'Datatype','int32');
    ncwrite(filename,'trinodes',Mesh.trinodes)
    % neighbours
    nccreate(filename,'nbe','Dimensions',{'nele' length(Mesh.uvnode) 'three' 3},'Datatype','int32');
    ncwrite(filename,'nbe',Mesh.nbe)
    % depth at uvnode
    nccreate(filename,'depthUvnode','Dimensions',{'elems' length(Mesh.uvdepth)},'Datatype','single');
    ncwrite(filename,'depthUvnode',Mesh.uvdepth)
    % depth at nodexy
    nccreate(filename,'depthNodexy','Dimensions',{'node' length(Mesh.depth)},'Datatype','single');
    ncwrite(filename,'depthNodexy',Mesh.depth)
    % siglay
    if ~isfield(Mesh,'sigvec')
        Mesh.sigvec=Mesh.siglay;
    end
    nccreate(filename,'siglay','Dimensions',{'siglay' length(Mesh.sigvec)},'Datatype','single');
    %ncwrite(filename, 'siglay',Mesh.siglay(1:10)) % removed 08/08/2019
    if (switchSigvec == 1)
        ncwrite(filename, 'siglay',-Mesh.sigvec)
    else
        ncwrite(filename, 'siglay',Mesh.sigvec)
    end
    
    % allBoundary
    if isfield(Mesh,'ISONB')
        nccreate(filename,'boundaryNodesAll','Dimensions',{'bnode' length(find(Mesh.ISONB))},'Datatype','int32');
        ncwrite(filename,'boundaryNodesAll',find(Mesh.ISONB))
    elseif isfield(Mesh,'isonb')
        nccreate(filename,'boundaryNodesAll','Dimensions',{'bnode' length(find(Mesh.isonb))},'Datatype','int32');
        ncwrite(filename,'boundaryNodesAll',find(Mesh.isonb))
    elseif isfield(Mesh,'boundaryNodesAll')
        nccreate(filename,'boundaryNodesAll','Dimensions',{'bnode' length(find(Mesh.boundaryNodesAll))},'Datatype','int32');
        ncwrite(filename,'boundaryNodesAll',find(Mesh.boundaryNodesAll))
    else
        try
            nccreate(filename,'boundaryNodesAll','Dimensions',{'bnode' length(Mesh.obc)+length(Mesh.closedBnodes)},'Datatype','int32');
            ncwrite(filename,'boundaryNodesAll',[Mesh.obc';Mesh.closedBnodes'])
        catch
            warning('Mesh does not contain field ISONB, obc or closedBnodes; did not write boundaryNodesAll');
        end
    end
    
    % openBoundary
    
    if ~isempty(obcFile)
        load(obcFile)
        nccreate(filename,'boundaryNodesOpen','Dimensions',{'obcnode' length(O.obc_nodes)},'Datatype','int32');
        ncwrite(filename,'boundaryNodesOpen',O.obc_nodes)
    else
        try
            nccreate(filename,'boundaryNodesOpen','Dimensions',{'obcnode' length(Mesh.obc)},'Datatype','int32');
            ncwrite(filename,'boundaryNodesOpen',Mesh.obc')
        catch
            warning('OBC data not available in Mesh file; did not write boundaryNodesOpen');
        end 
    end
    
    
    
    
    % % MeshType
    % nccreate(filename,'MeshType','Datatype','char');
    % ncwrite(filename,'MeshType','triangular');

    % dlmwrite('WestCOMS_boundaryAll.dat',find(Mesh.isonb))
    % dlmwrite('WestCOMS_boundaryOpen.dat',O.obc_nodes)
end