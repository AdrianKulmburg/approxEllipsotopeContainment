function [averageRuntime, stdError] = re_computeInputs_zonoLinSys(termRegObject, simulations, timeStepSize, N)

% Loading algorithm parameters
Nsimulation = length(simulations);

% Pre-processing
steps_indices = cell([1 Nsimulation]);
for i=1:Nsimulation
    current_index_list = [];
    for j=1:N
        current_index = find(simulations{i}.t >= (j-1)*timeStepSize, 1);
        current_index_list = [current_index_list current_index];
    end
    steps_indices{i} = current_index_list;
end



inputSets = termRegObject.inputSets;

% Suppress solver output
persistent options_linprog
if isempty(options_linprog)
    options_linprog = optimoptions('linprog', 'Display', 'none');
end

% And now we can start the show:
runtimes = zeros([Nsimulation 1]);
for i=1:Nsimulation
    currentTic = tic;
    y0_no_input = simulations{i}.y(:,1);
    % determine zonotope factors for input correction
    alpha = getInputScalings(termRegObject,y0_no_input,options_linprog);
        
    %y0_iter = simulations{i}.y(:,1);
    %u_iter = simulations{i}.u(:,1);
    
    for j=1:N
        % compute control input
        c_u = center(inputSets{j});
        G_u = generators(inputSets{j});
        G_u = G_u(:,1:length(alpha));
        
        y0_iter = simulations{i}.y(:,steps_indices{i}(j));
        u_iter = termRegObject.uEq + termRegObject.K*(y0_iter - termRegObject.xEq) + c_u + G_u*alpha;
    end
    runtimes(i) = toc(currentTic);
end


averageRuntime = sum(runtimes)/Nsimulation;
stdError = std(runtimes)/sqrt(Nsimulation);


end


function scalings = getInputScalings(obj,x0,options_linprog)
% compute scaling factors for the input zonotope
    
    G = obj.set.generators;
    c = obj.set.center;
    
    %%
    % See also cora/contSet/@zonotope/zonotopeNorm
    % Retrieve dimensions of the generator matrix
    n = size(G, 1);
    m = size(G, 2);
    
    % Set up objective and constraints of the linear program
    f = [1;sparse(m,1)];
    
    Aeq = [sparse(n, 1) G];
    beq = x0 - c;
    
    Aineq1 = [-ones([m 1]) speye(m)];
    Aineq2 = [-ones([m 1]) -speye(m)];
    
    Aineq = [Aineq1; Aineq2];
    bineq = sparse(2*m, 1);
    
    % Solve the linear program
    [preimage, fval, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, [], [], options_linprog);
    
    tol = 1e-5;
    
    % If the problem is not feasible, y has been chosen incorrectly
    if exitflag == -2 || fval > 1 + tol
        error('The point x does not lie in the corresponding reachable set.')
        % In case anything else went wrong, throw out an error
    elseif exitflag ~= 1
        throw(CORAerror('CORA:solverIssue'));
    end
    
    scalings = preimage(2:(m+1));
    
end