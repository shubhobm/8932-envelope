function [WX,W] = fusing(Y,X,u,hK)

for(h = 3:hK)
    i=h-2;
    parameters = read_inputs('lad','nslices',h);
    %[W,f] = lad(Y,X,1,'cont',parameters);
    Ys = slices(Y,parameters.nslices);
    data_parameters = setdatapars_v2(Ys,X,parameters.nslices);
    param(i)=data_parameters;
end
    
    Fhandle = F(@F4FLAD,param);
    dFhandle = dF(@dF4FLAD,param);

    
    p = cols(X); Wn = eye(p);
    fp = Fhandle(Wn);
    if u == p,
        disp('WARNING: the subspace you are looking for has the same dimension as the original feature space')
        fn = fp;
    else
        %--- get initial estimate .................................................
        if isempty(parameters.initvalue)||ischar(parameters.initvalue)
            guess = get_initial_estimate(Ys,X,u,data_parameters,parameters);
            Wo = guess(Fhandle);
        else
            Wo = parameters.initvalue;
        end
    
        if ~isempty(parameters.sg),
            [fn Wn] = sg_min(Fhandle,dFhandle,Wo,parameters.sg{:},parameters.maxiter);
        else
            [fn Wn] = sg_min(Fhandle,dFhandle,Wo,'prcg','euclidean',{1:u},'quiet',parameters.maxiter);
        end
    Wn = orth(Wn);
    end
    W=Wn;
    WX=X*W;
