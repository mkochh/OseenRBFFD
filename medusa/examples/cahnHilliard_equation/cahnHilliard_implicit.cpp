#include <medusa/Medusa.hpp>
#include <Eigen/SparseLU>
#include <random>
#include <iostream>
#include "BiharmonicOp.hpp"

int main() {
    /// Physical constants
    double D = 1;  // diffusion coefficient
    double L = 0.1;  // square of transition length
    double border = 10;

    /// Solver settings
    double dx = 0.2;
    double dt = 1e-2;
    double tEnd = 3;
    int neighbourNumber = 100;
    int monomialPower = 5;
    constexpr int dim = 2;
    double periodicBandWidth = 6;  // width of periodic band in units of dx

    /// Initial conditions
    int randomSeed = 110;

    /// Create the square domain and fill it with nodes.
    typedef mm::Vec<double, dim> vec_t;
    mm::BoxShape<vec_t> box(-border, border);
    mm::DomainDiscretization<vec_t> domain = box.discretizeBoundaryWithStep(dx);
    mm::GeneralFill<vec_t> fill;
    fill.seed(randomSeed);
    domain.fill(fill, dx);

    /// Create additional boundary nodes to satisfy periodic BC.
    // remove edges at + border and add edges at - border to interior
    Eigen::Matrix2d borders;  // col1: - border, col2: + border, rows for dimensions
    borders << domain.shape().bbox().first, domain.shape().bbox().second;
    for (auto idx : domain.boundary()) {
        const auto& pos = domain.pos(idx);
        if (pos[0] != borders(0, 1) && pos[1] != borders(1, 1)) {
            domain.changeToInterior(idx, pos, 1);
        }
    }
    domain.removeBoundaryNodes();

    // create new periodic boundary nodes for interior nodes that are closer to the edge than
    // periodicBandWidth * dx and construct a mapping between new and original nodes.
    mm::Range<int> periodicNodes;
    mm::Range<int> interiorMapping;
    mm::Vec2d normal(0, 0);  // not used in this case and can be chosen arbitrarily
    for (auto idx : domain.interior()) {
        const auto& pos = domain.pos(idx);
        // col1: distance to - border, col2: distance to + border, rows represent dimensions
        auto borderDistance = borders.colwise() - pos;
        const auto& isPeriodic = borderDistance.array().abs() < periodicBandWidth * dx;
        for (int i0 = 0; i0 < 2; i0++) {
            // check for periodicity in 1st dimension (x)
            if (isPeriodic(0, i0)) {
                auto periodicPos = pos;
                periodicPos(0) = borders(0, (i0 + 1) % 2) - borderDistance(0, i0);
                periodicNodes.push_back(domain.addBoundaryNode(periodicPos, -2, normal));
                interiorMapping.push_back(idx);
            }
            for (int i1 = 0; i1 < 2; i1++) {
                // check for periodicity in 2nd dimension (y)
                if (i0 == 0 && isPeriodic(1, i1)) {
                    auto periodicPos = pos;
                    periodicPos(1) = borders(1, (i1 + 1) % 2) - borderDistance(1, i1);
                    periodicNodes.push_back(domain.addBoundaryNode(periodicPos, -2, normal));
                    interiorMapping.push_back(idx);
                }
                // check for periodicity in both dimensions (corners)
                if (isPeriodic(0, i0) && isPeriodic(1, i1)) {
                    auto periodicPos = pos;
                    periodicPos(0) = borders(0, (i0 + 1) % 2) - borderDistance(0, i0);
                    periodicPos(1) = borders(1, (i1 + 1) % 2) - borderDistance(1, i1);
                    periodicNodes.push_back(domain.addBoundaryNode(periodicPos, -2, normal));
                    interiorMapping.push_back(idx);
                }
            }
        }
    }
    // copies for efficient iteration, domain.interior() is O(N) operation
    const auto& interiorNodes = domain.interior();
    const auto& allNodes = domain.all();
    int N = domain.size();
    int interiorSize = interiorNodes.size();
    int periodicSize = periodicNodes.size();

    /// Find support for the nodes and compute shapes.
    domain.findSupport(mm::FindClosest(neighbourNumber));

    mm::WLS<mm::Monomials<vec_t>, mm::GaussianWeight<vec_t>, mm::ScaleToClosest>
            approx(monomialPower);

    std::tuple<mm::Lap<dim>, Biharmonic<dim>> ops;
    mm::RaggedShapeStorage<vec_t, decltype(ops)> storage;
    storage.resize(domain.supportSizes());
    computeShapes(domain, approx, domain.all(),
                  std::tuple<mm::Lap<dim>, Biharmonic<dim>>(), &storage);

    /// Create system matrix and operators.
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    M.reserve(storage.supportSizes());
    Eigen::VectorXd rhs(N);
    rhs.setZero();

    auto expOp = storage.explicitOperators();
    auto impOp = storage.implicitOperators(M, rhs);

    /// Set initial concentration to random values.
    Eigen::VectorXd concentration(N);
    double initialSum = 0;
    std::mt19937 randomGenerator(randomSeed);
    std::uniform_real_distribution<double> realRandom(-1.0, 1.0);
    for (int i : interiorNodes) {
        concentration[i] = realRandom(randomGenerator);
        initialSum += concentration[i];
    }
    for (int i = 0; i < periodicSize; i++) {
        concentration[periodicNodes[i]] = concentration[interiorMapping[i]];
    }
    std::cout << "Initial average concentration: " << initialSum / interiorSize << std::endl;

    /// Create output file.
    int variableIdx = 0;
    mm::Range<std::string> matlabVariables = {"initial", "inter1", "inter2", "final"};
    mm::Range<double> printoutTimes = {0};
    std::ofstream out_file("cahnHilliard_visualization_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "border = " << border << ";" << std::endl;
    out_file << matlabVariables[variableIdx++] << " = " << concentration << ";" << std::endl;


    /// Set equations on internal and periodic nodes.
    for (int i : interiorNodes) {
        impOp.value(i) + D * dt * impOp.lap(i)
        + D * dt * L * impOp.apply<Biharmonic<dim>>(i) = rhs[i];
    }
    for (int i = 0; i < periodicSize; i++) {
        impOp.value(periodicNodes[i]) = rhs[interiorMapping[i]];
    }

    Eigen::SparseLU<decltype(M)> solver;
    solver.compute(M);

    /// Time stepping.
    double t = 0;
    Eigen::VectorXd cubedConcentration;
    while (t <= tEnd) {
        /// Prepare rhs.
        cubedConcentration = concentration.array() * concentration.array() * concentration.array();
        for (int i : allNodes) {
            concentration[i] = D * dt * expOp.lap(cubedConcentration, i) + concentration[i];
        }

        concentration = solver.solve(concentration);
        t += dt;

        /// Copy values from interior nodes to their periodic representations.
        for (int i = 0; i < periodicSize; i++) {
            concentration[periodicNodes[i]] = concentration[interiorMapping[i]];
        }

        if (std::fmod(t, tEnd / 3) < dt) {
            std::cout << "Printout at t = " << t
                      << ". Average concentration: "
                      << concentration(domain.interior()).sum() / interiorSize
                      << std::endl;
            out_file << matlabVariables[variableIdx++] << " = "
                     << concentration << ";" << std::endl;
            printoutTimes.push_back(t);
        }
    }
    out_file << "variables = " << matlabVariables << ";" << std::endl;
    out_file << "times = " << printoutTimes << ";" << std::endl;
    out_file.close();


    double finalSum = concentration(interiorNodes).sum();
    std::cout << std::endl
              << "Initial average concentration: " << initialSum / interiorSize << std::endl
              << "Final average concentration: " << finalSum / interiorSize << std::endl;

    return 0;
}
