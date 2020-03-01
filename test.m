function [R_r] = NR_Pipes()
    % parameters
    r = [1, 5, 1, 1, 5];
    % Initial conditions
    x = [1, 1, 1, 1, 1];

    for i = 1:100
        R = Compute_R(x, r, 2);
        R_r(i,:) = R;
        J = Compute_J(x, r, 2);
        dx = inv(J)*R';
        x = x - dx';
    end
    disp(x);
    semilogy(R_r.^2);
end

function [R] = Compute_R(x, r, n)
    R(1) = r(1)*x(1).^n + r(3)*x(3).^n - r(2)*x(2).^n;
    R(2) = r(5)*x(5).^n - r(4)*x(4).^n - r(3)*x(3).^n;
    R(3) = x(1) + x(2) - 10;
    R(4) = x(2) + x(3) - x(4);
    R(5) = x(1) - x(3) - x(5);
end

function [J] = Compute_J(x, r, n)
    J = [
        n*r(1)*x(1), -n*r(2)*x(2), n*r(3)*x(3), 0, 0;
        0, 0, -n*r(3)*x(3), -n*r(4)*x(4), n*r(5)*x(5);
        1, 1, 0, 0, 0;
        0, 1, 1, -1, 0;
        1, 0, -1, 0, -1
    ];
end
 
