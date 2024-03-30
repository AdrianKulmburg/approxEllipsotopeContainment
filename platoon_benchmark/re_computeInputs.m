function averageRuntime = re_computeInputs(termRegObject, simulations)

% Loading algorithm parameters
timeStepSize = termRegObject.optsInternal.Opts.timeStep;
N = termRegObject.optsInternal.Opts.N;
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

persistent options_quadprog
if isempty(options_quadprog)
    options_quadprog = optimoptions('quadprog', 'Display', 'none');
end

% Suppress solver output
persistent options_linprog
if isempty(options_linprog)
    options_linprog = optimoptions('linprog', 'Display', 'none');
end

% And now we can start the show:
totalRuntime = tic;
for i=1:Nsimulation
    y0_no_input = simulations{i}.y(:,1) - termRegObject.optsInternal.sys.D * simulations{i}.u(:,1);
    if strcmp(termRegObject.optsInternal.controlMethod, 'feedback')
        % for the feedback method, we need to compute the parametrization of
        % the initial point

        if strcmp(termRegObject.Opts.terminalRegionType, 'ellipsoid')
            if termRegObject.feedback_speedup.is_bijective
                % If the terminal region is just an ellipsoid, things are
                % really very easy
                preimage = termRegObject.feedback_speedup.multFactor * y0_no_input + termRegObject.feedback_speedup.constFactor;

                alpha = preimage(1:termRegObject.optsInternal.size_s);
                beta = preimage((termRegObject.optsInternal.size_s+1):end);
            else
                % Otherwise, we have to compute a slightly costly quadratic program
                V_ngen = size(generators(termRegObject.optsInternal.Param.V),2);

                H = 2.*blkdiag(speye(termRegObject.optsInternal.size_s), sparse(V_ngen,V_ngen));
                f = sparse(termRegObject.optsInternal.size_s+V_ngen,1);

                A = [sparse(V_ngen,termRegObject.optsInternal.size_s) speye(V_ngen); sparse(V_ngen, termRegObject.optsInternal.size_s) -speye(V_ngen)];
                b = ones([2*V_ngen 1]);

                Aeq = termRegObject.feedback_speedup.comb_G;
                beq = sparse(y0_no_input - termRegObject.feedback_speedup.comb_c);


                
                [preimage, fval, exitflag] = quadprog(H,f,A,b,Aeq,beq,[],[],options_quadprog);

                tol = 1e-5;

                % If the problem is not feasible, y has been chosen incorrectly
                if exitflag == -2 || fval > 1 + tol
                    error('The point x does not lie in the corresponding reachable set.')
                    % In case anything else went wrong, throw out an error
                elseif exitflag ~= 1
                    throw(CORAerror('CORA:solverIssue'));
                end

                alpha = preimage(1:termRegObject.optsInternal.size_s);
                beta = preimage((termRegObject.optsInternal.size_s+1):end);
            end

        else
            

            % See also cora/contSet/@zonotope/zonotopeNorm
            % Retrieve dimensions of the generator matrix
            n = size(termRegObject.feedback_speedup.G_comb, 1);
            m = size(termRegObject.feedback_speedup.G_comb, 2);

            % Set up objective and constraints of the linear program
            f = [1;sparse(m,1)];

            Aeq = [sparse(n, 1) termRegObject.feedback_speedup.G_comb];
            beq = simulations{i}.y(:,1) - termRegObject.optsInternal.sys.D*simulations{i}.u(:,1) - termRegObject.feedback_speedup.c_comb;

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

            alpha = preimage(2:(m - size(generators(termRegObject.optsInternal.Param.V),2)+1));
            beta = preimage((m - size(generators(termRegObject.optsInternal.Param.V),2)+2):end);
        end
    else
        alpha = [];
        beta = [];
    end
    %y0_iter = simulations{i}.y(:,1);
    u_iter = simulations{i}.u(:,1);
    for j=1:N
        y0_iter = simulations{i}.y(:,steps_indices{i}(j));
        u_iter = computeControlInput(termRegObject, y0_iter, u_iter, j, alpha, beta);
    end

end
totalRuntime = toc(totalRuntime);

averageRuntime = totalRuntime/Nsimulation;


end