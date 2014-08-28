function [WX,W] = fusing_disc(Y,X,u)


    parameters = read_inputs('lad');
    %[W,f] = lad(Y,X,1,'cont',parameters);
    % Ys = mapdata(Y);
    Ys=Y;
    parameters.nslices = max(Ys);
    data_parameters = setdatapars_v2(Ys,X,parameters.nslices);
    param(1)=data_parameters;

    for(i=1:length(Ys))
        if Ys(i)==3,
            Ys(i)=2;
        end
    end
    parameters.nslices = max(Ys);
    data_parameters = setdatapars_v2(Ys,X,parameters.nslices);
    param(2)=data_parameters;
    
    Ys=Y;
    for(i=1:length(Ys))
        if Ys(i)==3,
            Ys(i)=1;
        end
    end
    parameters.nslices = max(Ys);
    data_parameters = setdatapars_v2(Ys,X,parameters.nslices);
    param(3)=data_parameters;
    
    Ys=Y;
    for(i=1:length(Ys))
        if Ys(i)==3,
            Ys(i)=1;
        else
            Ys(i)=2;
        end
    end
    parameters.nslices = max(Ys);
    data_parameters = setdatapars_v2(Ys,X,parameters.nslices);
    param(4)=data_parameters;
    
    
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
