#include "steady_single_phase.hpp"


namespace CMF
{

bool validEntry(
    const Mesh& mesh,
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
    const Mesh& mesh,
    const BC& bc)
{
    for (size_t i = 1; i < mesh.size() - 1; ++i)
    {
        const auto& P = mesh[i];    // Center
        const auto& N = mesh[i+1];  // North
        const auto& S = mesh[i-1];  // South

        // Linear interpolation of viscosity at faces
        const real_t mu_n = (1.0 - P.width / N.width) * P.viscosity
                                 + P.width / N.width  * N.viscosity;
        const real_t mu_s = (1.0 - P.width / S.width) * P.viscosity
                                 + P.width / S.width  * S.viscosity;

        SM_ELEMENT_D(A, i, i + 1) = + mu_n * 2.0 / (P.width + N.width);
        SM_ELEMENT_D(A, i, i)     = - mu_n * 2.0 / (P.width + N.width)
                                    - mu_s * 2.0 / (P.width + S.width);
        SM_ELEMENT_D(A, i, i - 1) = + mu_s * 2.0 / (P.width + S.width);

        NV_Ith_S(b, i) = bc.global.second * P.width;
    }

    switch (bc.lowerWall.first)
    {
        case WallBC::Velocity:
        {
            const auto& P = mesh[0];
            const auto& N = mesh[1];

            const real_t mu_n = (1.0 - P.width / N.width) * P.viscosity
                                     + P.width / N.width  * N.viscosity;
            const real_t mu_s = P.viscosity; // TODO - Double check

            SM_ELEMENT_D(A, 0, 0) = + 2.0 * mu_n / (P.width + N.width)
                                    + 2.0 * mu_s / (P.width);
            SM_ELEMENT_D(A, 0, 1) = - 2.0 * mu_n / (P.width + N.width);
            NV_Ith_S(b, 0) = bc.global.second * P.width
                           + 2.0 * mu_s * bc.lowerWall.second / P.width;
            break;
        }

        case WallBC::VelocityGradient:
        {
            const auto& P = mesh[0];
            const auto& N = mesh[1];

            const real_t mu_n = (1.0 - P.width / N.width) * P.viscosity
                                     + P.width / N.width  * N.viscosity;
            const real_t mu_s = P.viscosity; // TODO - Double check

            SM_ELEMENT_D(A, 0, 0) = - 2.0 * mu_n / (P.width + N.width);
            SM_ELEMENT_D(A, 0, 1) = + 2.0 * mu_n / (P.width + N.width);
            NV_Ith_S(b, 0) = bc.global.second * P.width
                           + 2.0 * mu_s * bc.lowerWall.second / P.width;
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
            const size_t N = mesh.size() - 1;
            const auto& P = mesh[N];
            const auto& S = mesh[N - 1];

            const real_t mu_n = P.viscosity; // TODO - Double check
            const real_t mu_s = (1.0 - P.width / S.width) * P.viscosity
                                     + P.width / S.width  * S.viscosity;

            SM_ELEMENT_D(A, N, N)   = + 2.0 * mu_n / (P.width)
                                      + 2.0 * mu_s / (P.width + S.width);
            SM_ELEMENT_D(A, N, N-1) = - 2.0 * mu_s / (P.width + S.width);
            NV_Ith_S(b, N) = bc.global.second * P.width
                           + 2.0 * mu_n * bc.upperWall.second / P.width;
            break;
        }

        case WallBC::VelocityGradient:
        {
            const size_t N = mesh.size() - 1;
            const auto& P = mesh[N];
            const auto& S = mesh[N - 1];

            const real_t mu_n = P.viscosity; // TODO - Double check
            const real_t mu_s = (1.0 - P.width / S.width) * P.viscosity
                                     + P.width / S.width  * S.viscosity;

            SM_ELEMENT_D(A, N, N)   = - 2.0 * mu_n / (P.width + S.width);
            SM_ELEMENT_D(A, N, N-1) = + 2.0 * mu_s / (P.width + S.width);
            NV_Ith_S(b, N) = bc.global.second * P.width
                           - 2.0 * mu_n * bc.upperWall.second / P.width;
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


void steadyChannelFlow(
    Mesh& mesh,
    BC bc
)
{
    if (!validEntry(mesh, bc))
    {
        const char* msg = "Invalid input";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }

    const size_t N = mesh.size();
    LinearSolver solver(N);

    fillSystem(solver.Matrix(), solver.RHS(), mesh, bc);
    
    solver.solve();
    mesh.assignSolution(solver.getSolutionVector());
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

#if 0
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
#endif
} // namespace CMF
