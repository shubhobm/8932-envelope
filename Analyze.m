function [] = Analyze(X,Y,pemod,nam,typ,plot)

pcX = X*pemod.loadings;
s = size(pcX,1);

fprintf('Var explained by PCs and envelope gain ratios\n');
[pemod.pctvars pemod.envMod.ratio sqrt(length(pcX))*pemod.envMod.beta./pemod.envMod.asySE]

fprintf('Significant t-ratios\n');
sigpos = find(abs(pemod.tRatios)>1.96);
[nam(sigpos) num2cell([typ(sigpos)' pemod.tRatios(sigpos)])]
fprintf('%d among %d predictors fould significant.\n',length(sigpos),length(nam));

fprintf('Top 10 significant predictors:\n');
[a b] = sort(abs(pemod.tRatios));
for i = 1:9
    fprintf('%s, ',char(nam(b(i))));
end
fprintf('%s.\n', char(nam(b(10))));

if plot==1
    hold on
    grid on
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    for i=1:s
    if Y(i)==0    
        scatter3(pcX(i,1),pcX(i,2),pcX(i,3),'r.')
    else
        scatter3(pcX(i,1),pcX(i,2),pcX(i,3),'b.')
    end
    end
    hold off
end
% if pemod.u<3
% [a b] = meshgrid(-10:.1:10, -10:.1:10);
% c = pemod.envMod.Gamma(1,1).*a + pemod.envMod.Gamma(2,1).*b;
%     surf(a,b,c)
% end

