
function init = getInitialValues(t, p)
init = zeros(p.Netot+p.Nitot,1);
init(1:2*p.Ne,1) = -75 + 5*rand(2*p.Ne,1); %Vs, Vd %16/1/23
init(p.Netot+1:p.Netot+p.Ni,1) = -75 + 5*rand(p.Ni,1); %Vi %16/1/23
init(3*p.Ne+1:4*p.Ne,1) = p.Naeq; %Na+
init(2*p.Ne+1:3*p.Ne,:) = 1e-3; %Ca2+ concentration

%init = init';
end