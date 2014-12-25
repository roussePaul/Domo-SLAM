% This function displays the map landmark as black circles
%
function drawLandmarkMap(M)

if nargin < 1
    disp('You need to supply the file with landmarks to display (id x y)')
    return
end


N = size(M,3);
for i=1:N
    plot(M(:,1,i),M(:,2,i),'g');
end

