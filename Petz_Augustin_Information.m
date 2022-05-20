clear; close all;
% Change to the path where QETLAB and MFT toolboxes you downloaded
% If you have already added the toolboxes into your library directory, you
% don't need the following two lines.
addpath(genpath('C:\EXP\QETLAB-0.9'))
addpath(genpath('C:\EXP\mft\mftoolbox'))
%

% Random density matrix 
random_num = 1024; d = 10;
A = [];
for i = 1:random_num
    A(:,:,i) = RandomDensityMatrix(d);
end

p = rand(1, random_num); 
p = p./ sum(p); % probability of A

% Parameters 
max_iter_minimizer = 20;
max_iter = 10;
alpha = 0.5;

% for Armijo 
Armijo_params.alpha_mirror = 10; 
Armijo_params.r = 0.5;   
Armijo_params.tau = 0.5;

% for Polyak 1
Polyak1_params.delta_t = 2.5;
Polyak1_params.delta = 0.00001;
Polyak1_params.rho = 1.25; 
Polyak1_params.beta = 0.75; 

% for Polyak 2
Polyak2_params.delta_l = 2.5;
Polyak2_params.B = 2;

path = 0; tl = 1; l = 1; % Polyak 2 

sigma0 = eye(d)./sum(ones(1, d)); % initial point

% minimizer
sigma = sigma0; 
for i = 1 : max_iter_minimizer  
    X = ['Armijomin', num2str(i)];
    disp(X);
    func_grad = grad(sigma , alpha, p, A);
    sigma = Armijo(sigma, func_grad, alpha, p, A, Armijo_params);
end
f_min = func_val(sigma , alpha, p, A); % minimum value


% Algorithms:
% Armijo
seq_mirror = func_val(sigma0, alpha, p, A)- f_min;
seq_t_mirror = 0;
sigma = sigma0;

for i = 1 : max_iter
    X = ['Armijo', num2str(i)];
    disp(X);
    tic
    func_grad = grad(sigma, alpha, p, A);    
    [sigma] = Armijo(sigma, func_grad, alpha, p, A, Armijo_params);
    t_mirror = toc;
    seq_mirror= [seq_mirror max(0,real(func_val(sigma, alpha, p, A)- f_min))];
    seq_t_mirror = [seq_t_mirror, seq_t_mirror(i) + t_mirror]; 
end

% Polyak 1
seq_Polyak = func_val(sigma0, alpha, p, A)- f_min;
seq_t_Polyak = 0;

f_min_temp1 = func_val(sigma0, alpha, p, A);
sigma = sigma0;

for i = 1 : max_iter
    X = ['Polyak1', num2str(i)];
    disp(X);
    tic
    func_grad = grad(sigma , alpha, p, A);
    
    [sigma, Polyak1_params.delta_t] = Polyak1(sigma, func_grad, alpha, p, A, f_min_temp1, Polyak1_params);
    t_Polyak = toc;
    seq_Polyak = [seq_Polyak max(0,real(func_val(sigma, alpha, p, A)- f_min))];
    
    f_min_temp1 = min(f_min_temp1, func_val(sigma, alpha, p, A)); 
    
    seq_t_Polyak = [seq_t_Polyak, seq_t_Polyak(i) + t_Polyak];
end

% Polyak 2
seq_Polyak2 = func_val(sigma0, alpha, p, A)- f_min;
seq_t_Polyak2 = 0;
seq_Polyak_main2 = func_val(sigma0, alpha, p, A);
sigma = sigma0;

for i = 1 : max_iter
    X = ['Polyak2', num2str(i)];
    disp(X);
    tic
    func_grad = grad(sigma , alpha, p, A);
    
    [sigma, Polyak2_params.delta_l, path, tl, l] = Polyak2(sigma, func_grad, alpha, p, A, seq_Polyak_main2, i, path, tl, l, Polyak2_params);
    t_Polyak2 = toc;
    seq_Polyak2 = [seq_Polyak2 max(0,real(func_val(sigma, alpha, p, A)- f_min))];
    seq_Polyak_main2 = [seq_Polyak_main2 func_val(sigma, alpha, p, A)];
    
    seq_t_Polyak2 = [seq_t_Polyak2, seq_t_Polyak2(i) + t_Polyak2];
end

% Fixed-point iteration (Augustin)
seq_heuristic = func_val(sigma0, alpha, p, A)- f_min;
seq_t_heuristic = 0;
sigma = sigma0;

for i = 1:max_iter
    X = ['Heuristic', num2str(i)];
    disp(X);
    tic
    func_grad = grad(sigma , alpha, p, A);
    sigma = heuristic(sigma, func_grad);
    t_heuristic = toc;
    seq_heuristic = [seq_heuristic max(0,real(func_val(sigma, alpha, p, A)- f_min))];
    seq_t_heuristic = [seq_t_heuristic, seq_t_heuristic(i) + t_heuristic];
end

% do not consider numerical errors whose values are smaller than 10^(-5)
I = find(seq_mirror < 10^(-5), 1);
seq_mirror(:, I+1:end) = [];
I2 = find(seq_Polyak < 10^(-5), 1);
seq_Polyak(:, I2+1:end) = [];
I3 = find(seq_Polyak2 < 10^(-5), 1);
seq_Polyak2(:, I3+1:end) = [];
I4 = find(seq_heuristic < 10^(-5), 1);
seq_heuristic(:, I4+1:end) = [];

% f-f_min v.s iteration
figure(1)
lw = 'linewidth';
linewidth = 2;
subplot(2,1,1)
semilogy(real(seq_mirror), lw, linewidth)
hold on
semilogy(real(seq_Polyak), lw, linewidth)
semilogy(real(seq_Polyak2), lw, linewidth)
semilogy(real(seq_heuristic), lw, linewidth)
xlim([1 max([length(seq_mirror), length(seq_Polyak), length(seq_Polyak2) ,length(seq_heuristic)])])
ylim([10^(-5), max([max(real(seq_mirror)), max(real(seq_Polyak)), max(real(seq_Polyak2)), max(real(seq_heuristic))])])
yticks([10^(-5) 10^(-4) 10^(-2)])

grid on
xlabel('number of iterations')
ylabel('f(x_t)-f*')
legend("Armijo", "Polyak", "Polyak2", "Augustin")
hold off

% f-f_min v.s elapsed time
subplot(2,1,2)
lw = 'linewidth';
linewidth = 2;
semilogy(seq_t_mirror(:,1:length(seq_mirror)), seq_mirror, lw, linewidth)
hold on
semilogy(seq_t_Polyak(:,1:length(seq_Polyak)), seq_Polyak, lw, linewidth)
semilogy(seq_t_Polyak2(:,1:length(seq_Polyak2)), seq_Polyak2, lw, linewidth)
semilogy(seq_t_heuristic(:, 1:length(seq_heuristic)), seq_heuristic, lw, linewidth)
ylim([10^(-5), max([max(real(seq_mirror)), max(real(seq_Polyak)), max(real(seq_Polyak2)), max(real(seq_heuristic))])])
yticks([10^(-5) 10^(-4) 10^(-2)])

grid on
xlabel('Elapsed time (s)')
ylabel('f(x_t)-f*')
legend("Armijo", "Polyak", "Polyak2", "Augustin")
hold off

function [output,delta_t] = Polyak1(sigma, func_grad, alpha, p, A, f_min_temp1, Polyak1_params)
    
    delta = Polyak1_params.delta; rho = Polyak1_params.rho;
    beta = Polyak1_params.beta; delta_t = Polyak1_params.delta_t;
    
    f_t = func_val(sigma , alpha, p, A);
    f_min =  f_min_temp1;
    
    schattern_infi = abs(eigs(func_grad,1));
    
    eta =  (f_t - (f_min - delta_t))/(schattern_infi^2);
    
    x_alpha = expm(logm(sigma) - eta * func_grad);
    x_alpha = x_alpha / trace(x_alpha);
      
    if func_val(x_alpha , alpha, p, A) <= f_min - delta_t
        delta_t = rho * delta_t;
    else
        delta_t = max(beta * delta_t, delta);
    end
    output = x_alpha;

end

function [output, delta_l, path, tl, l] = Polyak2(sigma, func_grad, alpha, p, A, seq_Polyak_main, iter, path, tl, l, Polyak2_params)
    
    B = Polyak2_params.B; delta_l = Polyak2_params.delta_l;
    f_t = func_val(sigma , alpha, p, A);
    f_min_tl =  min(seq_Polyak_main(1:tl(l)));
    
    
    if f_t <= f_min_tl - delta_l/2
        tl = [tl, iter]; path = 0; l = l+1;
        f_min_tl =  min(seq_Polyak_main(1:tl(l)));
        
        schattern_infi = abs(eigs(func_grad,1));
        eta =  (f_t - (f_min_tl - delta_l))/(schattern_infi^2);
    
        x_alpha = expm(logm(sigma) - eta * func_grad);
        x_alpha = x_alpha / trace(x_alpha);
        
        path = path + eta * schattern_infi;
    else
        if path > B
            tl = [tl, iter]; path = 0; delta_l = delta_l/2;
            l = l+1; f_min_tl =  min(seq_Polyak_main(1:tl(l)));
            
            schattern_infi = abs(eigs(func_grad,1));    
            eta =  (f_t - (f_min_tl - delta_l))/(schattern_infi^2);
    
            x_alpha = expm(logm(sigma) - eta * func_grad);
            x_alpha = x_alpha / trace(x_alpha);
        
            path = path + eta * schattern_infi;
        else
            schattern_infi = abs(eigs(func_grad,1));    
            eta =  (f_t - (f_min_tl - delta_l))/(schattern_infi^2);
    
            x_alpha = expm(logm(sigma) - eta * func_grad);
            x_alpha = x_alpha / trace(x_alpha);
        
            path = path + eta * schattern_infi;
        end
    end
    output = x_alpha;

end


function [output] = Armijo(sigma, func_grad, alpha, p, A, Armijo_params)
    
    alpha_mirror = Armijo_params.alpha_mirror; r = Armijo_params.r; tau = Armijo_params.tau; 
    do = true;
    
    while do
        x_alpha = expm(logm(sigma) - alpha_mirror * func_grad);
        x_alpha = x_alpha / trace(x_alpha);
        
        do = tau * trace(func_grad * ((x_alpha - sigma)')) + func_val(sigma, alpha, p, A) <= func_val(x_alpha , alpha, p, A) - 1e-15;
        alpha_mirror = alpha_mirror * r;
    end
    output = x_alpha;
end

function [output] = heuristic(sigma, func_grad)
    temp1 = sigma^(1/2);
    output = -temp1 * func_grad * temp1;
end

% Cauculate the objecive function value
function output = func_val(sigma , alpha, p, A) 
    output = 0;
    sigma_t = sigma^ (1-alpha);
    for x = 1:length(p)
        temp = A(:,:,x);        
        temp = temp ^ (alpha);    
        output = output + p(x) * log(trace(temp * sigma_t));
    end
    output = (1/(alpha-1)) * output;
end

% Cauculate the gradient of the objecive function 
function output = grad(sigma , alpha, p, A) 
    output = 0 ;
    sigma1 = sigma^(1-alpha);
    for x = 1:length(p)
        temp = A(:,:,x); 
        temp = temp^(alpha);
        [~,Frechet] = powerm_pade_fre( sigma, 1-alpha, temp);
        output = output + p(x).* (Frechet) ./ (trace(temp * sigma1));      
    end
    output = output/(alpha-1); 
end


