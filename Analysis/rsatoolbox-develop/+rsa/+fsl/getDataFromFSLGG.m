function [betas,residuals] = getDataFromFSLGG(userOptions,subjectName)

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
    
    if 
    res=niftiread(fullfile(curfeatsPath,'stats','res4d.nii.gz'));
    res1=res(repmat(logical(mask),[1 1 1 size(res,4)]));
    residuals=cat(2,residuals,reshape(res1,[length(res1)/size(res,4) size(res,4)]));
    %ress=cat(2,ress,res(boolean(mask)));
end

residuals=residuals';
betas=betas';

end
