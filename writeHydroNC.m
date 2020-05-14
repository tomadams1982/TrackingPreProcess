function writeHydroNC(filename,Mesh,hydroDat,varargin)

    %filename='testFile.nc';
    
    ncid = netcdf.create(filename,'NOCLOBBER');

    dimid1 = netcdf.defDim(ncid,'nele',size(Mesh.uvnode,1));
    dimid2 = netcdf.defDim(ncid,'node',size(Mesh.nodexy,1));
    %dimid3a = netcdf.defDim(ncid,'one',1);
    dimid3b = netcdf.defDim(ncid,'two',2);
    dimid3c = netcdf.defDim(ncid,'three',3);
    dimid4 = netcdf.defDim(ncid,'nlevels',size(hydroDat.u,2));
    dimid5 = netcdf.defDim(ncid,'nhours',size(hydroDat.u,3));
    
    varid_mesh = netcdf.defVar(ncid,'fvcom_mesh' ,'NC_INT', [ ]       );
%     varid_uvnode = netcdf.defVar(ncid,'uvnode','NC_BYTE',[dimid1 dimid3b]);
%     varid_nodexy = netcdf.defVar(ncid,'nodexy','NC_BYTE',[dimid2 dimid3b]);
%     varid_trinodes = netcdf.defVar(ncid,'trinodes','NC_BYTE',[dimid1 dimid3c]);
%     varid_nbe = netcdf.defVar(ncid,'nbe','NC_BYTE',[dimid1 dimid3c]);
%     varid_u = netcdf.defVar(ncid,'u','NC_BYTE',[dimid1 dimid4 dimid5]);
%     varid_v = netcdf.defVar(ncid,'v','NC_BYTE',[dimid1 dimid4 dimid5]);
%     varid_lon = netcdf.defVar(ncid,'lon','NC_BYTE',[dimid2]);
%     varid_lat = netcdf.defVar(ncid,'lat','NC_BYTE',[dimid2]);
%     varid_lonc = netcdf.defVar(ncid,'lonc','NC_BYTE',[dimid1]);
%     varid_latc = netcdf.defVar(ncid,'latc','NC_BYTE',[dimid1]);

    netcdf.putAtt(ncid,varid_mesh,'cf_role'               ,'mesh_topology');
    netcdf.putAtt(ncid,varid_mesh,'topology_dimension'    ,'2');
    netcdf.putAtt(ncid,varid_mesh,'node_coordinates'      ,'lon lat');
    netcdf.putAtt(ncid,varid_mesh,'face_coordinates'      ,'lonc latc');
    netcdf.putAtt(ncid,varid_mesh,'face_node_connectivity','trinodes');

    netcdf.endDef(ncid);

    % uvnode
     nccreate(filename,'uvnode','Dimensions',{'nele' length(Mesh.uvnode) 'two' 2},'Datatype','single');
     ncwrite(filename,'uvnode',Mesh.uvnode)
    % nodexy
    nccreate(filename,'nodexy','Dimensions',{'node' length(Mesh.nodexy) 'two' 2},'Datatype','single');
    ncwrite(filename,'nodexy',Mesh.nodexy)
    % trinodes
    nccreate(filename,'trinodes','Dimensions',{'nele' length(Mesh.uvnode) 'three' 3},'Datatype','int32');
    ncwrite(filename,'trinodes',Mesh.trinodes)
    % neighbours
    nccreate(filename,'nbe','Dimensions',{'nele' length(Mesh.uvnode) 'three' 3},'Datatype','int32');
    ncwrite(filename,'nbe',Mesh.nbe)
    % u
    nccreate(filename,'u','Dimensions',{'nele' length(Mesh.uvnode) 'nlevels' size(hydroDat.u,2) 'nhours' size(hydroDat.u,3)},'Datatype','single');
    ncwrite(filename,'u',hydroDat.u)
    % v
    nccreate(filename,'v','Dimensions',{'nele' length(Mesh.uvnode) 'nlevels' size(hydroDat.u,2) 'nhours' size(hydroDat.u,3)},'Datatype','single');
    ncwrite(filename,'v',hydroDat.v)

    % lon
    nccreate(filename,'lon','Dimensions',{'node' length(Mesh.nodexy)},'Datatype','single');
    ncwrite(filename,'lon',Mesh.nodexy(:,1))
    % lat
    nccreate(filename,'lat','Dimensions',{'node' length(Mesh.nodexy)},'Datatype','single');
    ncwrite(filename,'lat',Mesh.nodexy(:,2))
    % lonc
    nccreate(filename,'lonc','Dimensions',{'nele' length(Mesh.uvnode)},'Datatype','single');
    ncwrite(filename,'lonc',Mesh.uvnode(:,1))
    % latc
    nccreate(filename,'latc','Dimensions',{'nele' length(Mesh.uvnode)},'Datatype','single');
    ncwrite(filename,'latc',Mesh.uvnode(:,2))
    


    
%     netcdf.putVar(ncid,varid_uvnode,Mesh.uvnode)
%     netcdf.putVar(ncid,varid_nodexy,Mesh.nodexy)
%     netcdf.putVar(ncid,varid_trinodes,Mesh.trinodes)
%     netcdf.putVar(ncid,varid_nbe,Mesh.nbe)
%     netcdf.putVar(ncid,varid_u,hydroDat.u)
%     netcdf.putVar(ncid,varid_v,hydroDat.v)

%     netcdf.putVar(ncid,varid_lon,Mesh.nodexy(:,1))
%     netcdf.putVar(ncid,varid_lat,Mesh.nodexy(:,2))
%     netcdf.putVar(ncid,varid_lonc,Mesh.uvnode(:,1))
%     netcdf.putVar(ncid,varid_latc,Mesh.uvnode(:,2))

    netcdf.close(ncid)
end