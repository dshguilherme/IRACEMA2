function geometry_obj = bilinear_coons(C0u,C0v,C1u,C1v)
assert(any(C0u.eval_point(0) == C0v.eval_point(0)), "Boundaries do not match");
assert(any(C0u.eval_point(1) == C1v.eval_point(0)), "Boundaries do not match");
assert(any(C1u.eval_point(0) == C0v.eval_point(1)), "Boundaries do not match");
assert(any(C1u.eval_point(1) == C1v.eval_point(1)), "Boundaries do not match");
S_00 = C0u.eval_point(0);
S_01 = C0v.eval_point(1);
S_11 = C1u.eval_point(1);
S_10 = C1v.eval_point(0);
Tx = @(u,v) [1 u]*[S_00(1) S_01(1); S_10(1) S_11(1)]*[1; v];
Ty = @(u,v) [1 u]*[S_00(2) S_01(2); S_10(2) S_11(2)]*[1; v];
Tz = @(u,v) [1 u]*[S_00(3) S_01(3); S_10(3) S_11(3)]*[1; v];
T = @(u,v) [Tx(u,v) Ty(u,v) Tz(u,v)];

T_00 = T(0,0);
T_01 = T(0,1);
T_11 = T(1,1);
T_10 = T(1,0);

R1 = bs_ruled_surface(C0u,C1u);
R2 = bs_ruled_surface(C0v,C1v);

%% Equaling polynomial degree and knots
    pu = max(R1.p(1),R2.p(1));
    if pu-R1.p(1) > 0
        R1.degree_elevate(1,pu-R1.p(1));
    elseif pu-R2.p(1) >0
        R2.degree_elevate(1,pu-R2.p(1));
    end

    pv = max(R1.p(2),R2.p(2));
    if pv-R1.p(2) > 0
        R2.degree_elevate(1,pv-R1.p(2));
    elseif pv-R2.p(2) >0
        R2.degree_elevate(1,pv-R2.p(1));
    end

    U1 = R1.knots{1};
    U2 = R2.knots{1};

    idx = ~ismember(U1,U2,'legacy');
    knots_to_add = U1(idx);
    if any(knots_to_add)
        R2.knot_refine(knots_to_add,1)
    end

    idx = ~ismember(U2,U1,'legacy');
    knots_to_add = U2(idx);
    if any(knots_to_add)
        R1.knot_refine(knots_to_add,1)
    end

    V1 = R1.knots{2};
    V2 = R2.knots{2};

    idx = ~ismember(V1,V2,'legacy');
    knots_to_add = V1(idx);
    if any(knots_to_add)
        R2.knot_refine(knots_to_add,2)
    end

    idx = ~ismember(V2,V1,'legacy');
    knots_to_add = V2(idx);
    if any(knots_to_add)
        R1.knot_refine(knots_to_add,2)
    end
%%
    

end