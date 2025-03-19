function S = PauliI(N,i,a)

S = 1;
switch a
    case 1
        for j = 1:N
            if j ~= i
                S = kron(S,eye(2));
            else
                S = kron(S,[0,1;1,0]);
            end
        end
    case 2
        for j = 1:N
            if j ~= i
                S = kron(S,eye(2));
            else
                S = kron(S,1i * [0,-1;1,0]);
            end
        end
    case 3
        for j = 1:N
            if j ~= i
                S = kron(S,eye(2));
            else
                S = kron(S,[1,0;0,-1]);
            end
        end
        
end

end