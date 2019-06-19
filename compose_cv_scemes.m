function SCHEMES = compose_cv_scemes(CVPATH)
    animate = [ones(50,1);zeros(50,1)];
    nitems = numel(animate);
    if isempty(CVPATH)
        nschemes = 10;
        nfolds = 10;
        SCHEMES = zeros(nitems, nschemes);
        for iScheme = 1:nschemes
            c = cvpartition(animate,'KFold', nfolds);
            for iFold = 1:nfolds
                SCHEMES(:,iScheme) = SCHEMES(:,iScheme) + (test(c, iFold) * iFold);
            end
        end
    else
        if exist(CVPATH, 'file')
            load(CVPATH, 'CV');
        else
            cvfile = CVPATH;
            CVPATH = fullfile(CV_DIR, cvfile);
            load(CVPATH, 'CV');
        end
        SCHEMES = CV;
    end
end