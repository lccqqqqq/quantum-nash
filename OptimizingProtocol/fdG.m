function dG = fdG(N,i,a)
% to create i[sigma,Z], the perturbed longitudinal field
dG = 1i*(PauliI(N,i,a)*PauliI(N,i,3)-PauliI(N,i,3)*PauliI(N,i,a));
end