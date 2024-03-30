function averageRuntime = re_computeInputs_zonoLinSys(termRegObject, simulations, timeStepSize, N)

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

% And now we can start the show:
totalRuntime = tic;
for i=1:Nsimulation
    y0_no_input = simulations{i}.y(:,1);
    % determine zonotope factors for input correction
    alpha = getInputScalings(termRegObject,y0_no_input);
        
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

end
totalRuntime = toc(totalRuntime);

averageRuntime = totalRuntime/Nsimulation;


end


function scalings = getInputScalings(obj,x0)
% compute scaling factors for the input zonotope

    % initialization
    center = obj.set.center;
    generators = obj.set.generators;
    [nx, nGen] = size(generators);
    
    % setup optimization problem
    state = sdpvar(nx, 1, 'full');
    scalings = sdpvar(nGen, 1, 'full');
    constraings = [state == center + generators * scalings];
    cost = norm(scalings, 'inf');
    % Suppress solver output
    persistent options_linprog
    if isempty(options_linprog)
        options_linprog = sdpsettings('verbose',0,'allownonconvex',0, ...
                                                    'solver','linprog');
    end
    opt = optimizer(constraings, cost, options_linprog, state, scalings);
    
    % compute input scaling by solving the optimization problem
    scalings = opt(x0);
end