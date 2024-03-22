clear all;



dimension = 20;
lambda_value = 10;

%{
number_of_constraints = 1;
for i = 1:number_of_constraints
    a = randn(dimension,dimension);
    a = a + a';
    [V,D] = eig(a);
    a = a/sum(abs(diag(D)));

    witness{i} = a;
end
%}

witness{1} = readmatrix('witness/witness_1_dim20.txt');

%witness{1} = readmatrix('witness/witness_1.txt');

x_axes = [];
upper_value = [];
lower_value = [];
for i = 100:100
grid = grid_function(1/8,0.0001,0.01,lambda_value);
a = inner_approx(witness, grid, dimension, lambda_value);

b = outer_approx(witness, grid, dimension, lambda_value);


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


function a = inner_approx(witness, grid, dimension, lambda_value)
    end_value = max(grid);
    [gamma, beta] = lists_gb_inner(grid);
    number_of_tau = length(grid);
    minimal_value = min(grid);
    
    cvx_solver Mosek
    cvx_begin sdp

    variable rho(dimension, dimension) hermitian semidefinite
    variable sigma(dimension, dimension) hermitian semidefinite
    variable X(dimension, dimension, number_of_tau) hermitian semidefinite
    
    for i = 1:number_of_tau
       X(:,:,i) >= beta(i)*rho + gamma(i)*sigma;
    end

    trace(rho) == 1;
    trace(sigma) == 1;
    rho <= lambda_value*sigma;
    minimal_value*sigma <= rho;
    for i = 1:size(witness,2)
        trace(witness{i}*rho)>=0.1;
        trace(witness{i}*sigma)<=-0.08;
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

function a = outer_approx(witness, grid, dimension, lambda_value)
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
       X(:,:,i) >= beta(i)*rho + gamma(i)*sigma;
    end

    trace(rho) == 1;
    trace(sigma) == 1;
    rho <= lambda_value*sigma;
    minimal_value*sigma <= rho;
    for i = 1:size(witness,2)
        trace(witness{i}*rho)>=0.1;
        trace(witness{i}*sigma)<=-0.08;
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

