function [alphan,fn,Wn] = new(Yaux,X,dim ,morph,parameters)
% [Wn,fn,fp] = new(Y,X,u,morph)
%=====================================================
u = dim(1); d=dim(2);

%----checking type of response ......................
if strcmpi(morph,'disc'),
    Y = mapdata(Yaux);
    parameters.nslices = max(Y);
else % morph = 'cont'
    Y = Yaux;
    parameters.nslices = length(Y);
end        

%--- get sample statistics ................................................
moreparameters = setdatapars_v2(Y,X,parameters.nslices);

if strcmpi(morph,'cont')
    [SIGMAfit,r] = get_fitted_cov(Y,X,parameters.fy);
else
    [SIGMAfit,r] = get_average_cov(X,moreparameters);
end
SIGMA = moreparameters.sigmag;
SIGMAres = SIGMA - SIGMAfit;
moreparameters.Afit = SIGMAfit;
moreparameters.B = SIGMAres;
moreparameters.u = u;
moreparameters.d = d;
moreparameters.r = r;



%--- get handle to objective function and derivative ......................
Fhandle = F(@F4new,moreparameters);
dFhandle = dF(@dF4new,Fhandle);

%--- get initial estimate .................................................
if strcmpi(morph,'cont')
    haux = 5;
    Ysliced = slices(Y,haux);
    aux_datapars = setdatapars_v2(Ysliced,X,haux);
    auxpars = parameters; auxpars.nslices=haux;
else
    Ysliced = Y;
    aux_datapars = moreparameters;
    auxpars = parameters; 
end
if isempty(parameters.initvalue)||ischar(parameters.initvalue)
    guess = get_initial_estimate(Ysliced,X,u,aux_datapars,auxpars);
    Wo = guess(Fhandle);
else
    Wo = parameters.initvalue;
end

%--- optimization .........................................................
p = cols(X); Wn = eye(p);
fp = Fhandle(Wn);
if u == p,
    disp('WARNING: the subspace you are looking for has the same dimension as the original feature space')
    fn = fp;
else
    if ~isempty(parameters.sg),
        [fn Wn] = sg_min(Fhandle,dFhandle,Wo,parameters.sg{:});
    else
        [fn Wn] = sg_min(Fhandle,dFhandle,Wo,'prcg','euclidean',{1:u},'quiet');
    end
end

aux=invsqrtm(Wn'*SIGMAres*Wn)*Wn'*SIGMA*Wn*invsqrtm(Wn'*SIGMAres*Wn);
[Vk,diagkD] = firsteigs(aux,d);
alphan=Wn*invsqrtm(Wn'*SIGMAres*Wn)*Vk;
alphan=orth(alphan);
