function W = pretpostavljena(x)

s = tf('s');
W =  tf(x(1)*[1/x(4) 1],conv([1/x(2) 1],[1/x(3) 1]));
%W =  (x(1)*exp(-x(4)*s))/((1/x(2)*s+1)*(1/x(3)*s+1));
%W =  tf(x(1),[1/x(2) 1]);
%W =  tf(x(1),conv([1/x(2) 1],[1/x(3) 1]));

