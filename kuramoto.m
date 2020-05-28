function dy=kuramoto(t,y)

nn=10; kk=5;

for i=1:nn
    w(i) = -1 + 2*(i-1)./(nn-1);
end

dy=zeros(nn,1);

% Flow equation

for i=1:nn
    
	dy(i) = w(i);
    
    
	for     j=1:nn
            
                dy(i) = dy(i) + (kk/nn)*sin(y(j)-y(i));
                
                disp((kk/nn)*sin(y(j)-y(i)))
    end
end

end