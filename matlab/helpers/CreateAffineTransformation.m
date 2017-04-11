function A = CreateAffineTransformation(trans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tx =  trans(1);
ty =  trans(2);
r2 =  trans(3);
sx =  trans(4);
sy =  trans(5);
r1 =  trans(6);


% syms tx ty r2 sx sy r1
% syms R1 T2 S A

R1 = [ cos(r1) , -sin(r1) ; sin(r1) , cos(r1) ];
R2 = [ cos(r2) , -sin(r2) ; sin(r2) , cos(r2) ];
S =  [ sx      , 0        ; 0       , sy      ];
A =  [0 0 tx ; 0 0 ty ; 0 0 1];
A(1:2,1:2) = R1*S*R2;


% A =     [ sx*cos(r1)*cos(r2) - sy*sin(r1)*sin(r2), - sx*cos(r1)*sin(r2) - sy*cos(r2)*sin(r1), tx]
%         [ sx*cos(r2)*sin(r1) + sy*cos(r1)*sin(r2),   sy*cos(r1)*cos(r2) - sx*sin(r1)*sin(r2), ty]
%         [                                       0,                                         0,  1]

% disp('done');


