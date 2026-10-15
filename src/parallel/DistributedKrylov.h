/**
 * @author Gabriel Henrique Silva
 * @file DistributedKrylov.h
 * @brief Row-distributed Jacobi-preconditioned BiCGSTAB (MPI).
 *
 * Complements the device-partitioned assembly: after the Allreduce fusion,
 * every rank holds the IDENTICAL Jacobian and RHS (row-replicated). This
 * solver distributes the Krylov work itself:
 *
 *  - each rank owns a contiguous ROW BLOCK [rowBegin, rowEnd);
 *  - sparse matvecs are evaluated per-block from the replicated matrix and
 *    the replicated input vector, then fused with ONE Allgatherv;
 *  - dot products / norms are partial sums over the owned block, fused
 *    with one scalar Allreduce;
 *  - all vectors (x, r, p, v, s, t) stay replicated and every rank executes
 *    the same scalar sequence, so the solution is bit-identical on all
 *    ranks with no solution broadcast needed.
 *
 * Cost per iteration: 2 × Allgatherv(n) + ~4 scalar Allreduces, vs. the
 * O(nnz) redundant matvec of the replicated solve. Worth it for large n;
 * small systems keep the replicated path (see Distribution thresholds).
 */

#ifndef DISTRIBUTED_KRYLOV_H
#define DISTRIBUTED_KRYLOV_H

#include <vector>

namespace e_XSpice {

class DaeMatrix;
class MPIManager;
struct SolverResult;

SolverResult distributedBiCGSTAB(const DaeMatrix& A,
                                 const std::vector<double>& b,
                                 std::vector<double>& x,
                                 int rowBegin, int rowEnd,
                                 MPIManager& mpi,
                                 int maxIterations,
                                 double tolerance);

} // namespace e_XSpice

#endif // DISTRIBUTED_KRYLOV_H
