function [hydroDat] = interpolateHydro(Mesh,data)
% INTERPOLATEHYDRO - interpolate a structure representing hydrodynamic
% model output from a ROMS model onto a triangular Mesh structure ("FVCOM"
% format)

    tic
    hydroDat.u = zeros(size(Mesh.uvnode,1),size(data.u,3),size(data.u,4));
    hydroDat.v = zeros(size(Mesh.uvnode,1),size(data.u,3),size(data.u,4));

    for dep=1:size(data.u,3)
        fprintf('depth %d : ',dep)
        for tt=1:size(data.u,4)
            fprintf('%d ',tt)
            Fu = scatteredInterpolant(reshape(data.lon_u,[],1),reshape(data.lat_u,[],1),reshape(squeeze(data.u(:,:,dep,tt)),[],1));
            hydroDat.u(:,dep,tt) = Fu(Mesh.uvnode(:,1),Mesh.uvnode(:,2));
            Fv = scatteredInterpolant(reshape(data.lon_v,[],1),reshape(data.lat_v,[],1),reshape(squeeze(data.v(:,:,dep,tt)),[],1));
            hydroDat.v(:,dep,tt) = Fv(Mesh.uvnode(:,1),Mesh.uvnode(:,2));
            
        end
        fprintf('\n')
    end
    
    hydroDat.u(isnan(hydroDat.u))=0;
    hydroDat.v(isnan(hydroDat.v))=0;
    
    toc

end