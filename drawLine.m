function drawLine(r,t)
x = linspace(-100,100,2);
y = (r-cos(t)*x)/sin(t);
plot(x,y,'c');
