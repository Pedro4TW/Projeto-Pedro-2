function res = consumption_equivalent(param,graph,c,L,welfare)
if param.mobility==1
    I=find(L>0);
    Uc=(c(I(1))/param.alpha)^param.alpha.*(param.Hj(I(1))/L(I(1))/(1-param.alpha))^(1-param.alpha);
    res=(welfare/Uc)^(1/param.alpha);
else
    if param.mobility==0
	omegaj=param.omegaj;
    elseif param.mobility==0.5
	omegaj=param.omegar(graph.region);
    end

    uj=param.u(c,param.Hj./L);
    uj(L==0)=0;
    Uc=sum(omegaj.*L.*uj);
    res=(welfare/Uc)^(1/(param.alpha*(1-param.rho)));
end
