# UnconstrainedOptimization_SimplexMethod

Non-linear Programming

Simplex Method\\
Chapter 6. (p.301) Non-Linear Programming II. Unconstrained Optimization
6.7 Simplex Method

Engineering Optimization - Theory and Practice, 4th Edition

Singiresu S. Rao

 The geometric figure formed by a set of (n+1) points in an n-dimensional
   space is called a simplex
 There three (3) parts of simplex:
   a). Reflection
           Xr = (1 + a)X0 - aXh
           Where:
               Xh := max f(Xi),    i = 1:n+1
               X0 := 1/n sum(Xi),  i = 1:n+1 & i =! h
   b). Expansion
           if f(Xr) < f(Xl), where Xl := min f(Xi),    i = 1:n+1
               Xe = gXr + (1 - g)X0
               if f(Xe) < f(Xl);   Xh = Xe (successful)
               if f(Xe) > f(Xl);   Xh = Xr (NOT successful)
           End if
   c). Contraction
           if f(Xr) > f(Xi),  i = 1:n+1 & i =! h
               if f(Xr) < f(Xh),
                   Xh = Xr
                   Xc = bXh + (1 - b)X0,
               elseif f(Xr) > f(Xh)
                   Xc = bXh + (1 - b)X0,
               End if
               if f(Xc) < min[f(Xh),f(Xr)],    Xh = Xc
                   otherwise, all Xi = (Xi + Xl)/2
               End if
           End if
