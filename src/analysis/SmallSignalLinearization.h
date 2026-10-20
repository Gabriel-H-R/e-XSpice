/**
 * @author Gabriel Henrique Silva
 * @file SmallSignalLinearization.h
 * @brief Shared preparation for the frequency-domain analyses (.AC, .NOISE):
 *        operating point -> linearization -> ground elimination ->
 *        pivot-safe ordering.
 *
 * Both frequency analyses need EXACTLY the same preparation, and it is the
 * part with all the subtlety (why the ground row is dropped, why the row
 * permutation is legal, why the pattern can be reused for every frequency).
 * Keeping it in one place is what makes "the noise analysis sees the same
 * system as the AC analysis" a structural fact instead of a coincidence.
 *
 * Sequence
 * --------
 *  1. DC operating point (SPICE/Xyce always compute one first), one DaeSystem
 *     and one nonlinear solve.
 *  2. refreshLinearization() so (G, C) belong to the CONVERGED point.
 *  3. Ground elimination: the DAE allocates a row/column for the ground node
 *     (with a tiny gmin on its diagonal) to keep the system square. For the
 *     frequency system ground is not an unknown at all, so its structurally
 *     empty row and the -g·x(gnd) entries of its column disappear — keeping
 *     them leaves a zero row that a no-pivot factorization can only survive
 *     by perturbing the diagonal, and the perturbation leaks into every other
 *     unknown (measured: 5e-7 error on a node an ideal source pins exactly).
 *  4. StaticRowPermutation: a ROW permutation that gives every column a usable
 *     diagonal (MNA rows like a node wired only to sources/inductors have
 *     none) with no fill-in, decided ONCE from the first frequency and frozen
 *     for the sweep.
 *
 * What a row permutation does NOT touch: x[j] belongs to COLUMN j, so the
 * unknown↔row mapping (acRow/fullRow) stays on the original indexing and
 * reporting is unaffected; only right-hand sides have to follow the rows,
 * which is what permuteRhs() does.
 */

#ifndef SMALL_SIGNAL_LINEARIZATION_H
#define SMALL_SIGNAL_LINEARIZATION_H

#include "../core/Circuit.h"
#include "../solver/ComplexSparseLU.h"

#include <memory>
#include <string>
#include <vector>

namespace e_XSpice {

class DaeSystem;

/**
 * @brief Operating point + (G, C) snapshot + ordering, ready for a sweep.
 */
class SmallSignalLinearization {
public:
    struct Options {
        int    max_iterations = 200;
        double rel_tol        = 1e-3;
        double abs_tol        = 1e-12;
        bool   verbose        = false;
        /** Solves the OP on the distributed device/row layout (MPI). */
        bool   use_distributed_layout = true;
        /**
         * @brief Frequency used to DECIDE the pivot-safe ordering (rad/s).
         *
         * The ordering only asks "does this column have a usable pivot?",
         * and the answer must not change over the sweep. For a frequency
         * sweep, any positive reference works; for the DC sensitivity the
         * matrix is G alone, so the reference must be 0 — otherwise an
         * inductor's branch column would look pivotable (j·ω·L on its
         * diagonal) when the DC matrix has no diagonal there at all.
         */
        double reference_omega = 6.28318530717958647692;   // 2*pi (1 Hz)
    };

    /**
     * @brief Runs the operating point and prepares the frequency system.
     * @return false on non-convergence or an unusable circuit; see message().
     */
    bool prepare(Circuit& circuit, const Options& opt);
    bool ok() const { return ok_; }
    const std::string& message() const { return message_; }

    /** Size of the ground-eliminated frequency system (the number of unknowns). */
    int size() const { return nn_; }
    /** Size of the DAE system (one more: the ground row/column). */
    int fullSize() const { return n_full_; }
    /** Rows 0..numNodeRows()-1 of the DAE are node voltages. */
    int numNodeRows() const { return num_node_rows_; }
    /** DAE row of the eliminated ground node (-1 when there is none). */
    int groundRow() const { return gnd_row_; }

    /** Frequency-system row of a DAE row (-1 for the ground row). */
    int acRow(int full_row) const {
        return (full_row >= 0 && full_row < n_full_)
                   ? ac_of_row_[static_cast<size_t>(full_row)]
                   : -1;
    }
    /** DAE row behind frequency-system row k. */
    int fullRow(int ac_row) const {
        return (ac_row >= 0 && ac_row < nn_)
                   ? row_of_ac_[static_cast<size_t>(ac_row)]
                   : -1;
    }

    // ---- the matrix, already filtered and permuted ----
    const std::vector<int>&    rowPointers() const { return row_ptr_; }
    const std::vector<int>&    colIndices()  const { return col_idx_; }
    const std::vector<double>& gValues()     const { return g_vals_; }
    const std::vector<double>& cValues()     const { return c_vals_; }

    const StaticRowPermutation& permutation() const { return perm_; }

    /**
     * @brief Maps a DAE-space vector onto the permuted frequency system.
     *
     * Used for right-hand sides, which live on the EQUATIONS (rows) and
     * therefore follow the permutation. Row k of the result is the component
     * of full row row_of_ac[perm.row_of_pos[k]].
     *
     * This is the correct transformation for a FORWARD solve (A'·x = b',
     * with A' = P·A): the right-hand side follows the rows.
     */
    std::vector<double> permuteRhs(const std::vector<double>& full) const;

    /**
     * @brief Drops the ground row (if any) KEEPING THE ORIGINAL ROW ORDER.
     *
     * This is the correct transformation for a TRANSPOSED solve
     * (A'^T·x = b): the transpose "undoes" the row permutation, so the
     * right-hand side must stay in the original indexing, and the SOLUTION
     * comes back in the permuted indexing (read it at acRow()).
     *
     * Getting this wrong is silent: the adjoint vector is then A^{-T}·P·c
     * instead of A^{-T}·c, which is right whenever the pivot-safe ordering
     * happens to be the identity and wrong for every circuit that needs a
     * row swap (a node wired only to sources/inductors, for instance).
     */
    std::vector<double> dropGroundRow(const std::vector<double>& full) const;

    /**
     * @brief Maps a permuted frequency-system vector back to DAE space.
     *
     * The inverse of permuteRhs(): the solution of the frequency system
     * lives in the permuted indexing, but the devices read their node
     * voltages in the DAE indexing, so the AC solution has to be mapped
     * back before it is handed to a device.
     */
    std::vector<double> unpermuteSolution(const std::vector<double>& perm) const;

    /**
     * @brief Maps a permuted frequency-system vector back to DAE space,
     *        writing into a caller-provided buffer (no allocation).
     */
    void unpermuteSolution(const std::vector<double>& perm,
                            std::vector<double>& out) const;

    /** Newton iterations spent on the operating point. */
    int newtonIterations() const { return newton_iterations_; }

    /** The converged operating point (valid while this object lives). */
    DaeSystem& system() { return *system_; }
    const DaeSystem& system() const { return *system_; }

private:
    bool ok_ = false;
    std::string message_;
    int n_full_ = 0;
    int num_node_rows_ = 0;
    int gnd_row_ = -1;
    int nn_ = 0;
    int newton_iterations_ = 0;

    std::vector<int> ac_of_row_;   // DAE row -> frequency row
    std::vector<int> row_of_ac_;   // frequency row -> DAE row

    // CSR of A = G + j·w·C: ground dropped AND rows permuted.
    std::vector<int>    row_ptr_;
    std::vector<int>    col_idx_;
    std::vector<double> g_vals_;
    std::vector<double> c_vals_;

    StaticRowPermutation perm_;

    // Owned (not RCP): the DAE system is not a RefCounted type in this
    // codebase, and it must outlive this object because the operating point
    // it holds is what the noise sources are evaluated at.
    std::unique_ptr<DaeSystem> system_;
    Circuit* circuit_ = nullptr;
};

} // namespace e_XSpice

#endif // SMALL_SIGNAL_LINEARIZATION_H