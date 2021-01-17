function [ cost ] = costFunction( del,maxJump,alpha,param )
%COSTFUNCTION calculate the cost between the nodes. The farther the nodes
%are, the bigger the cost will be.
%   MAXJUMP is the smallest distance from which the cost is applied (ie >0)
%   ALPHA is the parameter of the sigmoid. the biggest this number is, the
%   smoothes the sigmoid will be. The sigmoid is aimed at including non
%   linearity (ex: if the distance between the nodes is doubled, the
%   penality is not doubled)
%   param : method of the calculus. 
            %'non linear', add a non linear term. The penality increase with the square of the distance 
            %'linear' the penality is proportional with the distance multiplied with a sigmoid
            % to foster small distances. 
switch param
    case 'nonLinear'
        alphaNL = 1;
        sigmOffset =0.5;
        % linear part
        linear = (heaviside(del-maxJump) .* abs((del-maxJump))) .* (sigmf(del,[alpha,maxJump])-sigmOffset);
        % non linear part
        nonLinear = alphaNL.*(sigmf(del,[0.01,maxJump])-sigmOffset).*heaviside(del-maxJump).*abs(del-maxJump).^2;
        cost = linear + nonLinear;
        
    case 'linear'
        sigmOffset =0.5;
        % linear part
        linear = (heaviside(del-maxJump) .* abs((del-maxJump))) .* (sigmf(del,[alpha,maxJump])-sigmOffset);
        % non linear part
        cost = linear;
        
end

