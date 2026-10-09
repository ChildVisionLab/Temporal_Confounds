function [betas,residuals, nVols] = getDataFromFSL(userOptions)

subjectName = userOptions.subjectNames{:};
%load mask
curmaskPath=rsa.util.replaceWildcards(userOptions.maskPath, ...
        '[[subjectName]]',subjectName, ...
        '[[maskName]]',userOptions.maskNames{1});
mask=niftiread(curmaskPath);

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
        beta=niftiread(fullfile(curfeatsPath,'stats',['cope' int2str(c{1}) '.nii.gz']));
        betas=cat(2,betas,beta(logical(mask)));
    end
    %106   106    48   406
    res=niftiread(fullfile(curfeatsPath,'stats','res4d.nii.gz'));
    nVols(contains(userOptions.run_names, r)) = size(res,4);
    % 459998           1
    res1=res(repmat(logical(mask),[1 1 1 size(res,4)]));
    % 1133         406
    residuals=cat(2,residuals,reshape(res1,[length(res1)/size(res,4) size(res,4)]));
    %ress=cat(2,ress,res(boolean(mask)));
end

residuals=residuals';
betas=betas';

end
