function geometry_obj = coons(C0u,C0v,C1u,C1v,D0u,D1u,D0v,D1v)
% Calculate S_ij
assert(any(C0u.eval_point(0) == C0v.eval_point(0)), "Boundaries do not match");
assert(any(C0u.eval_point(1) == C1v.eval_point(0)), "Boundaries do not match");
assert(any(C1u.eval_point(0) == C0v.eval_point(1)), "Boundaries do not match");
assert(any(C1u.eval_point(1) == C1v.eval_point(1)), "Boundaries do not match");
S_00 = C0u.eval_point(0);
S_01 = C0v.eval_point(1);
S_11 = C1u.eval_point(1);
S_10 = C1v.eval_point(0);
% S_00 -> S_01 -> S_11 -> S_10

D0u = C0v;
D1u = C1v;
D0v = C0u;
D1v = C1u;

% Compute ctrl points of S1(u,v)
Pk0 = C0u.points;
Pk1 = C1u.points;
Qk0 = D0u.points;
Qk1 = D1u.points;

P_S1 = cell(length(Pk0),4);
P_S1(:,[1 2]) = Pk0(1);
P_S1(:,[end-1 end]) = Pk1(1);

    der = Qk0{1};
    der = der(:,[1:3]).*der(:,4);
    der = [der 0];
    der2 = Qk1{1};
    der2 = der2(:,[1:3]).*der2(:,4);
    der2 = [der2 0];
    for s=1:length(Pk0)
        P_S1{s,2} = P_S1{s,2} +(1/3)*der;
        P_S1{s,end-1} = P_S1{s,end-1} -(1/3)*der2;
    end

S1 = Geometry(2,{C0u.knots{1},[0 0 0 0 1 1 1 1]}, P_S1, [C0u.p(1) 3]);
    
% Compute ctrl points of S2(u,v)
Pl0 = C0v.points;
Pl1 = C1v.points;
Ql0 = D0v.points;
Ql1 = D1v.points;

P_S2 = cell(4,length(Pl0));
P_S2([1 2],:) = Pl0(1);
P_S2([end-1 end],:) = Pl1(1);

    der = Ql0{1};
    der = der(:,[1:3]).*der(:,4);
    der = [der 0];
    der2 = Ql1{1};
    der2 = der2(:,[1:3]).*der2(:,4);
    der2 = [der2 0];
    for s=1:length(Pl0)
        P_S2{2,s} = P_S2{2,s} +(1/3)*der;
        P_S2{end-1,s} = P_S2{end-1,s} -(1/3)*der2;
    end
    
S2 = Geometry(2,{[0 0 0 0 1 1 1 1], C0v.knots{1},}, P_S2, [3 C0v.p(1)]);

% Compute ctrl points of T(u,v)

P_T00 = cell(2,2);
P_T00{1,1} = [S_00 1];
P_T00{2,1} = [((1/3)*D0v.eval_point(0)+S_00) 1];
P_T00{1,2} = [(1/3)*D0u.eval_point(0)+S_00 1];
P_T00{2,2} = P_T00{2,1} +P_T00{1,2} -P_T00{1,1};

P_T01 = cell(2,2);
P_T01{1,1} = [S_01 1];
P_T01{2,1} = [(1/3)*D0v.eval_point(1)+S_01 1];
P_T01{1,2} = [(1/3)*D0u.eval_point(1)+S_01 1];
P_T01{2,2} = P_T01{2,1} +P_T01{1,2} -P_T01{1,1};

P_T11 = cell(2,2);
P_T11{1,1} = [S_11 1];
P_T11{2,1} = [(1/3)*D1v.eval_point(1)+S_11 1];
P_T11{1,2} = [(1/3)*D1u.eval_point(1)+S_11 1];
P_T11{2,2} = P_T11{2,1} +P_T11{1,2} -P_T11{1,1};

P_T10 = cell(2,2);
P_T10{1,1} = [S_10 1];
P_T10{2,1} = [(1/3)*D1v.eval_point(0)+S_10 1];
P_T10{1,2} = [(1/3)*D1u.eval_point(0)+S_10 1];
P_T10{2,2} = P_T10{2,1} +P_T10{1,2} -P_T10{1,1};
% S_00 -> S_01 -> S_11 -> S_10
P_T = cell(4,4);
P_T(1:2,1:2) = P_T00;
P_T(1:2,3:4) = P_T01;
P_T(3:4,1:2) = P_T10;
P_T(3:4,3:4) = P_T11;

T = Geometry(2,{[0 0 0 0 1 1 1 1], [0 0 0 0 1 1 1 1]}, P_T, [3 3]);

for i=1:2
pu = max(3,max(S1.p(i),S2.p(i)));
    if pu-S1.p(i) > 0
        S1.degree_elevate(pu-S1.p(i),i);
    elseif pu-S2.p(i) >0
        S2.degree_elevate(pu-S2.p(i),i);
    elseif pu-T.p(i) >0
        T.degree_elevate(pu-T.p(i),i);
    end
end
 
for i=1:2
    U1 = S1.knots{i};
    U2 = S2.knots{i};
    U3 = T.knots{i};
  
    idx = ~ismember(U1,U2,'legacy');
    knots_to_add = U1(idx);
    if any(knots_to_add)
        S2.knot_refine(knots_to_add,i)
    end

    idx = ~ismember(U2,U1,'legacy');
    knots_to_add = U2(idx);
    if any(knots_to_add)
        S1.knot_refine(knots_to_add,i)
    end
    
    idx = ~ismember(U1,U3,'legacy');
    knots_to_add = U1(idx);
    if any(knots_to_add)
        T.knot_refine(knots_to_add,i);
    end
    
    idx = ~ismember(U2,U3,'legacy');
    knots_to_add = U2(idx);
    if any(knots_to_add)
        T.knot_refine(knots_to_add,i);
    end
end
    
    
points = cell(size(T.points));
[i_points j_points] = size(points);
for i=1:i_points
    for j=1:j_points
        points{i,j} = S1.points{i,j} +S2.points{i,j} -T.points{i,j};
    end
end
geometry_obj = Geometry(2,T.knots,points,T.p);
end