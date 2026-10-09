function [betas,residuals, nVols] = getSurfDataFromFSL(userOptions, hemi)

subjectName = userOptions.subjectNames{:};

%load copes
betas=[];
residuals=[];
for r = userOptions.run_names
    curfeatsPath=rsa.util.replaceWildcards(userOptions.featsPath, ...
        '[[subjectName]]',subjectName, ...
        '[[featPrefix]]',userOptions.featsPrefix, ...
        '[[runName]]',r{1}, ...
        '[[featSuffix]]',userOptions.featsSuffix);
    for c = userOptions.copes 
        beta=fs_load_mgh(fullfile(curfeatsPath,'stats',['cope' int2str(c{1}) '_surf_' hemi '.mgh']));
        betas=cat(2,betas,beta);
    end
    %   131725           1           1         406
    res=fs_load_mgh(fullfile(curfeatsPath,'stats',['res4d_surf_', hemi, '.mgh']));
    res = squeeze(res);
    nVols(contains(userOptions.run_names, r)) = size(res,2);
    residuals=cat(2,residuals,res);
end
residuals=residuals';
betas=betas';
end
