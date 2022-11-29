%% Simplex Method
% Chapter 6. (p.301) Non-Linear Programming II. Unconstrained Optimization
% 6.7 Simplex Method
% Engineering Optimization - Theory and Practice, 4th Edition
% Singiresu S. Rao

% The geometric figure formed by a set of (n+1) points in an n-dimensional
%   space is called a simplex
% There three (3) parts of simplex:
%   a). Reflection
%           Xr = (1 + a)X0 - aXh
%           Where:
%               Xh := max f(Xi),    i = 1:n+1
%               X0 := 1/n sum(Xi),  i = 1:n+1 & i =! h
%   b). Expansion
%           if f(Xr) < f(Xl), where Xl := min f(Xi),    i = 1:n+1
%               Xe = gXr + (1 - g)X0
%               if f(Xe) < f(Xl);   Xh = Xe (successful)
%               if f(Xe) > f(Xl);   Xh = Xr (NOT successful)
%           End if
%   c). Contraction
%           if f(Xr) > f(Xi),  i = 1:n+1 & i =! h
%               if f(Xr) < f(Xh),
%                   Xh = Xr
%                   Xc = bXh + (1 - b)X0,
%               elseif f(Xr) > f(Xh)
%                   Xc = bXh + (1 - b)X0,
%               End if
%               if f(Xc) < min[f(Xh),f(Xr)],    Xh = Xc
%                   otherwise, all Xi = (Xi + Xl)/2
%               End if
%           End if

clear all; clc

%% Initialization
syms X Y lambda
var = 2;                    % Number of variables on objective function
n = 2;                      % The dimension
eps = 0.01;                 % Some prescribed small quantity (\epsilon)
a = 1.0;                    % Alpha
b = 0.5;                    % Beta
c = 2.0;                    % Gamma

x1 = [4;4];                 % Starting point X1
x2 = [5;4];                 % Starting point X2
x3 = [4;5];                 % Starting point X3

% Generate objective function
f = symfun(X - Y + 2*X^2 + 2*X*Y + Y^2, [X Y]);     % 1 Method
% f = inline(X - Y + 2*X^2 + 2*X*Y + Y^2)           % 2 Method
% f = @(X,Y) X - Y + 2*X^2 + 2*X*Y + Y^2            % 3 Method

Q = inf;                    % Initial standard deviation

% Constructing the set of (n+1) points in n-dimensional space
X = [x1,x2,x3];
F1 = f(X(1,1),X(2,1));
F2 = f(X(1,2),X(2,2));
F3 = f(X(1,3),X(2,3));
F = [F1,F2,F3];

%% Looping Simplex
while Q > eps

    % Finding vertex Xh from the maximum among F1,F2,F3
    [Fh, idMax] = max(F);                   Xh = X(:,idMax);
    % Finding the minimum Xl among F1,F2,F3
    [Fl, idMin] = min(F);                   Xl = X(:,idMin);

    % REFLECTION
    X_max = X;
    X_max(:,idMax) = [];
    X0 = (1/n)*sum(X_max')';                F0 = f(X0(1),X0(2));
    F0 = double(F0);
    Xr = (1+a)*X0 - a*Xh;                   Fr = f(Xr(1),Xr(2));
    Fr = double(Fr);

    % EXPANSION
    if Fr < Fl
        Xe = c*Xr + (1-c)*X0;               Fe = f(Xe(1),Xe(2));
        Fe = double(Fe);
        if Fe < Fl
            Xh = Xe;
        else
            Xh = Xr;
        end
    end

    % CONTRACTION
    F_temp = F;
    F_temp(idMax) = [];
    if Fr > F_temp
        if Fr < Fh
            Xh = Xr;
        end
        Xc = b*Xh + (1-b)*X0;               Fc = f(Xc(1),Xc(2));
        Fc = double(Fc);

        if Fc < min([Fh,Fr])
            Xh = Xc;
        else
            for i = 1:n+1
                % Update X if all failed
                X(i) = (X(i) + Xl)/2;
                % Update the functions of vertices based on new-Xi
                F(i) = f(X(1,i),X(2,i));
            end
        end
    end

    % Update the vertices X
    X(:,idMax) = Xh;
    Fh = f(Xh(1),Xh(2));
    % Update the functions of vertices based on new-X
    F(idMax) = Fh;

    % Computing the standard deviation
    Sum = (F(1)-F0).^2 + (F(2)-F0).^2 + (F(3)-F0).^2;
    Q = double(sqrt(Sum/(n+1)))

    % Plottng the Simplex
    axis square
    Z = linspace(-7.5,7.5);             Y = linspace(0,7.5);
    [A,B] = meshgrid(Z,Y);
    f_fig = f(A,B);
    levels = 10:10:350;
    figure(1), contour(Z,Y,f_fig,levels,'LineWidth',.5), colorbar
    hold on;
    Sx = [X X(:,1)];
    plot(Sx(1,:),Sx(2,:),'LineWidth',2);
end

% This code is dedicated for educational porposes
% The writer do encourage to users to try/take a look/check what's going on
%   for each step so that they can completely understand
% The complete algorithm and example could be found inside the book