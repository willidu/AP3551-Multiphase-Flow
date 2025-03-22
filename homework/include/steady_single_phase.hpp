#pragma once

#include "core.hpp"

namespace CMF
{
    
// ****************************************************************************
//                                    TASK 1
// ****************************************************************************
// Develop a routine to compute the velocity in a steady fully-developed 
// single-phase channel flow with variable viscosity, using the finite-volume 
// method. Assume that the viscosity profile is a prescribed input to the 
// routine.
//
// Consider three possible wall boundary conditions: 
//   (i)   prescribed velocity at the wall, 
//   (ii)  prescribed velocity-gradient at the wall, and 
//   (iii) prescribed relation between the shear-stress at the wall and the 
//         velocity field. 
//
// Consider two possible global boundary conditions: 
//   (i)  prescribed pressure-gradient and 
//   (ii) prescribed flow-rate.
// ****************************************************************************

enum class WallBC   { Velocity, VelocityGradient, ShearStress };
enum class GlobalBC { PressureGradient, FlowRate };

struct BC
{
    std::pair<WallBC, real_t> lowerWall;
    std::pair<WallBC, real_t> upperWall;
    std::pair<GlobalBC, real_t> global;
};

struct Mesh
{
    struct Node
    {
        real_t y;         // Node position
        real_t width;     // Cell width
        real_t u;         // Velocity
        real_t viscosity; // Effective (dynamic) viscosity

        Node() = default;
        Node(real_t y, real_t width, real_t u, real_t viscosity)
            : y(y), width(width), u(u), viscosity(viscosity) {}
    };

    Mesh() = default;
    Mesh(size_t N, real_t H, real_t viscosity)
    {
        nodes.reserve(N);
        real_t dy = H / N;
        for (size_t i = 0; i < N; ++i)
            nodes.emplace_back((i + 0.5) * dy, dy, 0.0, viscosity);
    }
    Mesh(std::vector<Node> nodes) : nodes(nodes) {}

    std::vector<Node> nodes;

    Node operator[](size_t i) const
    {
        return nodes.at(i);
    }

    size_t size() const noexcept { return nodes.size(); }

    void assignSolution(const std::vector<real_t>& u)
    {
        assert(u.size() == nodes.size() && "Invalid solution vector size");
        for (size_t i = 0; i < nodes.size(); ++i)
            nodes.at(i).u = u.at(i);
    }

    std::pair<std::vector<real_t>, std::vector<real_t>> getSolution() const
    {
        std::vector<real_t> y, u;
        y.reserve(nodes.size());
        u.reserve(nodes.size());
        for (const auto& node : nodes)
        {
            y.push_back(node.y);
            u.push_back(node.u);
        }
        return {y, u};
    }

    void applyViscosityProfile(std::function<real_t(real_t)> profile)
    {
        for (auto& node : nodes)
            node.viscosity = profile(node.y);
    }
};

/**
 * Find the velocity field in a steady fully-developed single-phase channel flow
 * with variable viscosity. The viscosity profile is a prescribed input to the
 * routine.
 */
void steadyChannelFlow(Mesh& mesh, BC bc);

# if 0
// ****************************************************************************
//                                    TASK 2
// ****************************************************************************
// Develop a routine to compute the eddy-viscosity Prandtl mixing-length model,
// with prescribed functions for the mixing-length (try different functions).
// Incorporate this routine into the channel flow code.
// ****************************************************************************

/**
 * @brief Mixing-length model for the eddy viscosity.
 * 
 * Find the velocity field in a steady fully-developed single-phase channel flow
 * with variable viscosity and mixing length. The viscosity profile is a
 * prescribed input to the routine. The mixing length is used to compute the
 * eddy viscosity.
 * 
 * @param mesh         Mesh points
 * @param viscosity    Molecular viscosity.
 * @param mixingLength The mixing length model.
 * @param bc           Boundary conditions.
 * @param maxIter      Maximum number of non-linear iterations.
 * @param tol          Stopping criteria for non-linear iterations.
 * 
 * @return Channel velocity field with eddy viscosity.
 */
std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    real_t viscosity,
    real_t density,
    std::function<real_t(real_t)> mixingLength,
    BC bc,
    size_t maxIter = 500,
    real_t tol = 1e-6
);


// ****************************************************************************
//                                    TASK 3
// ****************************************************************************
// Develop a routine to use a log-law wall-function as the prescribed relation
// between the shear-stress at the wall and the velocity-field. Consider the
// possibility of both smooth and rough walls. Incorporate this routine into
// the channel flow code.
// ****************************************************************************

std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    real_t viscosity,
    real_t density,
    std::function<real_t(real_t)> mixingLength,
    BC bc,
    real_t averageRoughness,
    std::function<real_t(real_t)> dampingFunction,
    size_t maxIter = 1000,
    real_t tol = 1e-6
);

#endif
} // namespace CMF
