clear all

addpath(genpath('../QETLAB-0.9'))

A1 = 1/sqrt(2) * [1,1;1,-1];
B1 = 1/sqrt(2) * [1,1;1i,-1i];

for i = 1:size(A1,1)
    for j = 1:size(B1,2)
        if abs(abs(A1(i,:)*B1(:,j))-1/sqrt(2)) > 0.01
            disp('error 2');
        end
    end
end

A2 = kron(A1,A1);
B2 = kron(B1,B1);

for i = 1:size(A2,1)
    for j = 1:size(B2,2)
        if abs(abs(A2(i,:)*B2(:,j))-1/sqrt(4)) > 0.01
            disp('error 4');
        end
    end
end

A3 = kron(A2,A1);
B3 = kron(B2,B1);

for i = 1:size(A3,1)
    for j = 1:size(B3,2)
        if abs(abs(A3(i,:)*B3(:,j))-1/sqrt(8)) > 0.01
            disp('error 8');
        end
    end
end

A4 = kron(A3,A1);
B4 = kron(B3,B1);

for i = 1:size(A4,1)
    for j = 1:size(B4,2)
        if abs(abs(A4(i,:)*B4(:,j))-1/sqrt(16)) > 0.01
            disp('error 16');
        end
    end
end






dimension = 64;
local_dimension = 8;
alpha = 0.1;
grid = grid_function(1/8,0.05,0.01,local_dimension);
register = measurement_statistic(alpha, local_dimension, A3, B3);

%a = inner_approx(register, A3, grid, dimension);
%b = outer_approx(register, A3, grid, dimension);










function a = inner_approx(register, basis, grid, dimension)
    end_value = max(grid);
    [gamma, beta] = lists_gb_inner(grid);
    number_of_tau = length(grid);
    
    
    cvx_solver Mosek
    cvx_begin sdp

    variable rho(dimension, dimension) hermitian semidefinite
    variable sigma(dimension, dimension) hermitian semidefinite
    variable X(dimension, dimension, number_of_tau) hermitian semidefinite
    
    for i = 1:number_of_tau
       X(:,:,i) >= beta(i)*rho + gamma(i)*pinching_alice(rho, basis);
    end

    trace(rho) == 1;
    

    for i = 1:size(register{2},2)
        trace(rho*register{2}{i}) == register{1}(i);
    end


  
    a = 0;
    for i = 1:number_of_tau
        a = a + trace(X(:,:,i));
    end

    % Define the objective function
    minimize(real(a));
    cvx_end 
    a = a + log(end_value) + 1 - end_value
    
end

function a = outer_approx(register, basis, grid, dimension)
    end_value = max(grid);
    minimal_value = min(grid);
    [gamma, beta] = lists_gb_outer(grid);
    number_of_tau = length(grid) - 1;

    cvx_precision high
    cvx_solver Mosek
    cvx_begin sdp

    variable rho(dimension, dimension) hermitian semidefinite
    variable sigma(dimension, dimension) hermitian semidefinite
    variable X(dimension, dimension, number_of_tau) hermitian semidefinite
    
    for i = 1:number_of_tau
       X(:,:,i) >= beta(i)*rho + gamma(i)*pinching_alice(rho, basis);
    end

    trace(rho) == 1;
    

    for i = 1:size(register{2},2)
        trace(rho*register{2}{i}) == register{1}(i);
    end


    a = 0;
    for i = 1:number_of_tau
        a = a + trace(X(:,:,i));
    end
   
   
    % Define the objective function
    minimize(real(a));
    cvx_end 
    a = a + log(end_value) + 1 - end_value
end














function grid = grid_function(c, epsilon, mu, lamb)
    grid = [mu];
    f = mu;
    i = 1;
    while f < lamb
        f = grid(i) + sqrt(grid(i)*epsilon/c);
        grid = [grid, f];
        i = i + 1;
    end
    grid(end) = lamb;
end





function [gamma, beta] = lists_gb_outer(grid)
    gamma = zeros(1, length(grid)-1);
    beta = zeros(1, length(grid)-1);
    for i = 1:length(grid)-1
        gamma(i) = grid(i+1) - grid(i);
        beta(i) = log(grid(i) / grid(i+1));
    end
end

function [gamma, beta] = lists_gb_inner(grid)
    gamma = zeros(1, length(grid));
    beta = zeros(1, length(grid));
    for i = 1:length(grid)
        if i == 1
            constante = (1 + grid(i) / (grid(i+1) - grid(i))) * log(grid(i+1) / grid(i)) - 1;
            gamma(i) = grid(i) * constante;
            beta(i) = -constante;
        elseif i == length(grid)
            constante = 1 - grid(i-1) / (grid(i) - grid(i-1)) * log(grid(i) / grid(i-1));
            gamma(i) = grid(i) * constante;
            beta(i) = -constante;
        else
            constante = (1 + grid(i) / (grid(i+1) - grid(i))) * log(grid(i+1) / grid(i)) - grid(i-1) / (grid(i) - grid(i-1)) * log(grid(i) / grid(i-1));
            gamma(i) = constante * grid(i);
            beta(i) = -constante;
        end
    end
end


function output_state = pinching_alice(state, basis)
    dimension = size(basis, 2);
    output_state = zeros(dimension^2, dimension^2);
    
    for i = 1:size(basis,2)
        element = basis(:,i);
        projector = element * element';
        kraus_operator = kron(projector, eye(dimension));
        output_state = output_state + kraus_operator * state * kraus_operator';
    end
end


function phi_matrix = phi(alpha, dimension)
    % Compute maximally mixed state
    maximally_mixed = (1/(dimension^2)) * eye(dimension^2);
    
    % Initialize phi_state
    phi_state = zeros(1, dimension^2);
    
    % Compute phi_state
    for i = 1:dimension
        a = zeros(1, dimension);
        a(i) = 1;
        phi_state = phi_state + kron(a, a);
    end
    
    % Normalize phi_state
    phi_state = (1/sqrt(dimension)) * phi_state;
    
    % Compute phi_matrix
    phi_matrix = (1-alpha) * (phi_state' * phi_state) + alpha * maximally_mixed;
end



function register = measurement_statistic(alpha, dimension, basis1, basis2)
    projections = {};
    for i = 1:size(basis1,2)
        projections{i} = basis1(:,i) * basis1(:,i)';
    end

    for j = 1:size(basis2,2)
        projections{j} = basis2(:,j) * basis2(:,j)';
    end
    
    probabilities = [];
    j = 1;
    for k = 1:size(projections,2)
        for l = 1:size(projections,2)
            operator = kron(projections{k}, projections{l});
            meas_operatoren{j} = operator;
            probabilities = [probabilities,trace(operator * phi(alpha, dimension))];
            j = j + 1;
        end
    end
    
    register{1} = real(probabilities);
    register{2} = meas_operatoren;
end

