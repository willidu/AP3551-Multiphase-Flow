#include "steady_single_phase.hpp"


namespace CMF
{

bool validEntry(
    const std::vector<real_t>& mesh,
    const BC& bc)
{
    if (mesh.size() < 3)
    {
        LOG_ERROR("Mesh must have at least 3 nodes");
        return false;
    }

    if (bc.global.first != GlobalBC::PressureGradient)
    {
        LOG_ERROR("Only PressureGradient BC is supported");
        return false;
    }

    if (bc.lowerWall.first == WallBC::VelocityGradient
        && bc.upperWall.first == WallBC::VelocityGradient)
    {
        LOG_ERROR("Both walls cannot have VelocityGradient BC");
        return false;
    }

    return true;
}


void fillSystem(
    SUNMatrix A,
    N_Vector b,
    const std::vector<real_t>& mesh,
    const std::function<real_t(real_t)>& viscosity,
    const BC& bc)
{
    // TODO: Check indexes (i, i+1/2, i-1/2 etc.)

    const size_t N = mesh.size();

    for (size_t i = 1; i < N - 1; ++i)
    {
        // Forward and backward spacing
        const real_t dy_f = mesh.at(i+1) - mesh.at(i);
        const real_t dy_b = mesh.at(i) - mesh.at(i-1);
        const real_t dy = (dy_f + dy_b) / 2.0;

        const real_t mu_ip = viscosity((mesh.at(i) + mesh.at(i+1)) / 2);
        const real_t mu_im = viscosity((mesh.at(i) + mesh.at(i-1)) / 2);

        SM_ELEMENT_D(A, i, i - 1) = + mu_im / (dy_b * dy);
        SM_ELEMENT_D(A, i, i)     = - mu_ip / (dy_f * dy)
                                    - mu_im / (dy_b * dy);
        SM_ELEMENT_D(A, i, i + 1) = + mu_ip / (dy_f * dy);

        NV_Ith_S(b, i) = bc.global.second;
    }

    switch (bc.lowerWall.first)
    {
        case WallBC::Velocity:
        {
            const real_t dy = mesh.at(1) - mesh.at(0);
            const real_t mu0 = viscosity(mesh.at(0) + dy / 2.0);
            SM_ELEMENT_D(A, 0, 0) = 3 * mu0 / (dy * dy);
            SM_ELEMENT_D(A, 0, 1) = -mu0 / (dy * dy);
            NV_Ith_S(b, 0) = bc.global.second 
                           + (2 * mu0 * bc.lowerWall.second / (dy * dy));
            break;
        }

        case WallBC::VelocityGradient:
        {
            const real_t dy = mesh.at(1) - mesh.at(0);
            const real_t mu0 = viscosity(mesh.at(0) + dy / 2.0);
            SM_ELEMENT_D(A, 0, 0) = 2 * mu0 / (dy * dy);
            SM_ELEMENT_D(A, 0, 1) = -2 * mu0 / (dy * dy);
            NV_Ith_S(b, 0) = bc.global.second
                           - bc.lowerWall.second / dy;
            break;
        }

        default:
        {
            const char* msg = "Unsupported lower wall BC";
            LOG_ERROR(msg);
            throw std::runtime_error(msg);
        }
    }

    switch (bc.upperWall.first)
    {
        case WallBC::Velocity:
        {
            const real_t dy = mesh.at(N-1) - mesh.at(N-2);
            const real_t muN = viscosity(mesh.at(N-1) - dy / 2.0);
            SM_ELEMENT_D(A, N - 1, N - 2) = -muN / (dy * dy);
            SM_ELEMENT_D(A, N - 1, N - 1) = 3 * muN / (dy * dy);
            NV_Ith_S(b, N - 1) = bc.global.second
                               + (2 * muN * bc.upperWall.second / (dy * dy));
            break;
        }

        case WallBC::VelocityGradient:
        {
            const real_t dy = mesh.at(N-1) - mesh.at(N-2);
            const real_t muN = viscosity(mesh.at(N-1) - dy / 2.0);
            SM_ELEMENT_D(A, N - 1, N - 2) = -2 * muN / (dy * dy);
            SM_ELEMENT_D(A, N - 1, N - 1) = 2 * muN / (dy * dy);
            NV_Ith_S(b, N - 1) = bc.global.second
                               - bc.upperWall.second / dy;
            break;
        }

        default:
        {
            const char* msg = "Unsupported upper wall BC";
            LOG_ERROR(msg);
            throw std::runtime_error(msg);
        }
    }
}


std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    std::function<real_t(real_t)> viscosity,
    BC bc
)
{
    assert(viscosity && "Invalid viscosity function");

    if (!validEntry(mesh, bc))
    {
        const char* msg = "Invalid input";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }

    const size_t N = mesh.size();
    LinearSolver solver(N);

    fillSystem(solver.Matrix(), solver.RHS(), mesh, viscosity, bc);
    
    solver.solve();
    return solver.getSolutionVector();
}


size_t findIndex(
    const std::vector<real_t>& mesh,
    real_t y
)
{
    const auto it = std::lower_bound(mesh.begin(), mesh.end(), y);

    if (it == mesh.end())
        return mesh.size() - 1;

    if (it == mesh.begin())
        return 0;

    return std::distance(mesh.begin(), it);
}


real_t velocityGradient(
    const std::vector<real_t>& mesh,
    const std::vector<real_t>& u,
    size_t i
)
{
    if (i == 0)
        return (u.at(1) - u.at(0)) / (mesh.at(1) - mesh.at(0));
    
    if (i == mesh.size() - 1)
        return (u.at(i) - u.at(i-1)) / (mesh.at(i) - mesh.at(i-1));

    return (u.at(i+1) - u.at(i-1)) / (mesh.at(i+1) - mesh.at(i-1));
}


std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    real_t viscosity_mol,
    real_t density,
    std::function<real_t(real_t)> mixingLength,
    BC bc,
    size_t maxIter,
    real_t tol
)
{
    assert(viscosity_mol > 0.0 && "Invalid molecular viscosity");
    assert(density > 0.0 && "Invalid density");

    if (!validEntry(mesh, bc))
    {
        const char* msg = "Invalid input";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }

    const size_t N = mesh.size();

    // We will fist assume that the eddy viscosity is zero
    std::function<real_t(real_t)> viscosity = [&](real_t y) -> real_t
    {
        return viscosity_mol;
    };

    std::vector<real_t> u;
    {
        LinearSolver solver(N);
        fillSystem(solver.Matrix(), solver.RHS(), mesh, viscosity, bc);
        solver.solve();
        u = solver.getSolutionVector();
    }

    // We will now iterate to find a solution with non-zero eddy viscosity
    for (size_t iter = 0; iter < maxIter; ++iter)
    {
        LinearSolver solver(N);

        viscosity = [&](real_t y) -> real_t
        {
            return viscosity_mol + density
                                 * std::pow(mixingLength(y), 2.0) 
                                 * std::abs(velocityGradient(mesh, u, findIndex(mesh, y)));
        };

        fillSystem(solver.Matrix(), solver.RHS(), mesh, viscosity, bc);
        solver.solve();
        std::vector<real_t> u_new = solver.getSolutionVector();

        // Convergence with L2 norm
        const real_t error = std::sqrt(std::inner_product(
            u.begin(), u.end(), u_new.begin(), 0.0, std::plus<real_t>(),
            [](real_t a, real_t b) { return std::pow(a - b, 2.0); }
        ));

        u = u_new;

        if (error < tol)
        {
            LOG_INFO("Converged after {0} iterations", iter);
            break;
        }

        if (iter == maxIter - 1)
        {
            LOG_WARN(
                "Did not converge after {0} iterations. Final residual: {1:.2e}",
                maxIter, error);
        }
    }

    return u;
}


std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    real_t viscosity_mol,
    real_t density,
    std::function<real_t(real_t)> mixingLength,
    BC bc,
    real_t averageRoughness,
    std::function<real_t(real_t)> dampingFunction,
    size_t maxIter,
    real_t tol
)
{
    assert(viscosity_mol > 0.0 && "Invalid molecular viscosity");
    assert(density > 0.0 && "Invalid density");
    assert(averageRoughness >= 0.0 && "Invalid average roughness");

    if (!validEntry(mesh, bc))
    {
        const char* msg = "Invalid input";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }

    const size_t N = mesh.size();

    // We will fist assume that the eddy viscosity is zero
    std::function<real_t(real_t)> viscosity = [&](real_t y) -> real_t
    {
        return viscosity_mol;
    };

    std::vector<real_t> u;
    {
        LinearSolver solver(N);
        fillSystem(solver.Matrix(), solver.RHS(), mesh, viscosity, bc);
        solver.solve();
        u = solver.getSolutionVector();
    }

    // We will now iterate to find a solution with non-zero eddy viscosity
    for (size_t iter = 0; iter < maxIter; ++iter)
    {
        LinearSolver solver(N);

        viscosity = [&](real_t y) -> real_t
        {
            auto yplus_utau = [&](real_t y) -> std::pair<real_t, real_t>
            {
                if (y < 0.5 * mesh.back()) // Bottom half
                {
                    const real_t dudy = velocityGradient(mesh, u, 0);
                    const real_t tau = viscosity_mol * std::abs(dudy);
                    const real_t u_tau = std::sqrt(tau / density);
                    const real_t yplus = y * u_tau / (viscosity_mol / density);
                    return {yplus, u_tau};
                }
                else // Top half
                {
                    const real_t dudy = velocityGradient(mesh, u, mesh.size()-1);
                    const real_t tau = viscosity_mol * std::abs(dudy);
                    const real_t u_tau = std::sqrt(tau / density);
                    const real_t yplus = (mesh.back() - y) * u_tau / (viscosity_mol / density);
                    return {yplus, u_tau};
                }
            };

            const auto [yplus, u_tau] = yplus_utau(y);
            // LOG_TRACE("y = {2:.3f}, yplus = {0:.1f}, u_tau = {1:.2f}", yplus, u_tau, y);

            if (yplus < 3.0)
            {
                // -> U+ = y+
                return viscosity_mol + density * 0.41 * y * yplus;
            }
            else if (yplus < 20.0)
            {
                // LOG_WARN("Buffer layer");
            }
            else if (yplus > 30.0 && yplus < 800.0)
            {
                // Log-law layer
                static constexpr real_t kappa = 0.41;
                static constexpr real_t beta = 5.2;

                return viscosity_mol + density * kappa * y * u_tau; // ?

                /*
                    H = 2*delta
                    tau = 2*nu*rho*abs(U[0])/dy
                    u_tau = np.sqrt(tau/rho)
                    y_plus = np.append(y[y<=delta]*u_tau/nu,(H-y[y>delta])*u_tau/nu)
                    fmu = (1-np.exp(-y_plus/A))**2
                */


            }

            return viscosity_mol + density
                                 * std::pow(mixingLength(y), 2.0) 
                                 * std::abs(velocityGradient(mesh, u, findIndex(mesh, y)));
                                //  * dampingFunction(yplus);
        };

        fillSystem(solver.Matrix(), solver.RHS(), mesh, viscosity, bc);
        solver.solve();
        std::vector<real_t> u_new = solver.getSolutionVector();

        // Convergence with L2 norm
        const real_t error = std::sqrt(std::inner_product(
            u.begin(), u.end(), u_new.begin(), 0.0, std::plus<real_t>(),
            [](real_t a, real_t b) { return std::pow(a - b, 2.0); }
        ));

        u = u_new;

        if (error < tol)
        {
            LOG_INFO("Converged after {0} iterations", iter);
            break;
        }

        if (iter == maxIter - 1)
        {
            LOG_WARN("Did not converge after {0} iterations", maxIter);
        }
    }

    return u;
}

} // namespace CMF
