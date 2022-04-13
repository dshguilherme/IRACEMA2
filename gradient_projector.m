function solution = gradient_projector(domain, fun)
asb = Poisson(1,"gauss",3,domain);
F1 = asb.project_grad(fun);
B1 = asb.L2_projector(3);
lm = asb.location_matrix;
id = asb.id_matrix;

% Dirichlet boundary conditions for projectors

boundaries = domain.extract_boundaries;

end