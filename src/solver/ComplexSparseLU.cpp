/**
 * @author Gabriel Henrique Silva
 * @file ComplexSparseLU.cpp
 * @brief Sparse complex LU: A = L·U, L unit diagonal, static pattern.
 *
 * Left-looking (column-oriented) Crout, no pivoting, natural column order:
 *
 *     column k:  x = A(:,k)
 *                for i < k:  u_ik = x[i]          (U(i,k), already known)
 *                           x -= u_ik · L(:,i)     (column i of L)
 *                L(k,k) = x[k]  and  L(i,k) = x[i]/L(k,k)  for i > k
 *
 * The update walks COLUMNS of L, which is why the combined L|U pattern is
 * kept column-ordered with the diagonal split; the forward substitution then
 * needs no division (unit L) and the backward one divides only by U(k,k).
 */

#include "ComplexSparseLU.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstdlib>

namespace e_XSpice {

// ===============================================================
// Static pivot-safe row permutation
// ===============================================================

StaticRowPermutation StaticRowPermutation::build(const std::vector<int>& rowPtr,
                                                 const std::vector<int>& colIdx,
                                                 const std::vector<double>& re,
                                                 const std::vector<double>& im,
                                                 const std::vector<int>& priority) {
    StaticRowPermutation perm;
    const int n = static_cast<int>(rowPtr.size()) - 1;
    perm.row_of_pos.assign(static_cast<size_t>(std::max(0, n)), -1);
    if (n <= 0) {
        perm.ok = true;
        return perm;
    }
    std::vector<char> taken(static_cast<size_t>(n), 0);

    // Magnitude of every entry, and the largest one of each column: a column
    // can carry its own pivot only if its diagonal is not negligible compared
    // to the couplings it actually has.
    const size_t nnz = colIdx.size();
    std::vector<double> mag(nnz, 0.0);
    std::vector<double> col_max(static_cast<size_t>(n), 0.0);
    for (size_t p = 0; p < nnz && p < re.size() && p < im.size(); ++p) {
        mag[p] = std::hypot(re[p], im[p]);
        const int c = colIdx[p];
        if (c >= 0 && c < n) col_max[static_cast<size_t>(c)] =
            std::max(col_max[static_cast<size_t>(c)], mag[p]);
    }
    // A diagonal below this fraction of the column maximum counts as "not
    // there" (1e-8 keeps legitimate small conductances, rejects the zeros and
    // the denormal-ish residues of a cancelled stamp).
    const double kUsableRatio = 1e-8;
    auto diagonalUsable = [&](int k) -> bool {
        for (int p = rowPtr[k]; p < rowPtr[k + 1]; ++p) {
            if (colIdx[static_cast<size_t>(p)] != k) continue;
            const double d = mag[static_cast<size_t>(p)];
            const double m = col_max[static_cast<size_t>(k)];
            return m > 0.0 && d >= kUsableRatio * m;
        }
        return false;   // no diagonal entry at all
    };

    // Pass 1: keep every usable natural diagonal (never worth trading away).
    for (int k = 0; k < n; ++k) {
        if (diagonalUsable(k) && !taken[static_cast<size_t>(k)]) {
            perm.row_of_pos[static_cast<size_t>(k)] = k;
            taken[static_cast<size_t>(k)] = 1;
        }
    }

    // Pass 2: columns that cannot carry their own pivot take the best free row
    // of their pattern: branch rows first (priority), then the strongest
    // coupling, with the index as a deterministic tie-break.
    const auto prio = [&](int r) -> int {
        return (r >= 0 && r < static_cast<int>(priority.size())) ? priority[static_cast<size_t>(r)]
                                                                 : 0;
    };
    for (int k = 0; k < n; ++k) {
        if (perm.row_of_pos[static_cast<size_t>(k)] >= 0) continue;
        ++perm.weak_diagonals;
        int best = -1;
        double best_mag = -1.0;
        for (int p = rowPtr[k]; p < rowPtr[k + 1]; ++p) {
            const int r = colIdx[static_cast<size_t>(p)];
            if (r < 0 || r >= n || taken[static_cast<size_t>(r)]) continue;
            const double m = mag[static_cast<size_t>(p)];
            if (best < 0 || prio(r) > prio(best) ||
                (prio(r) == prio(best) && (m > best_mag ||
                                           (m == best_mag && r < best)))) {
                best = r;
                best_mag = m;
            }
        }
        if (best < 0) {
            perm.ok = false;
            perm.reason = "no pivot-safe ordering exists for column " +
                          std::to_string(k) +
                          " (structurally singular matrix)";
            return perm;
        }
        perm.row_of_pos[static_cast<size_t>(k)] = best;
        taken[static_cast<size_t>(best)] = 1;
        if (best != k) ++perm.swaps;
    }
    perm.ok = true;
    return perm;
}

ComplexSparseLU::ComplexSparseLU() = default;

void ComplexSparseLU::setPattern(const std::vector<int>& rowPtr,
                                 const std::vector<int>& colIdx) {
    row_ptr_ = rowPtr;
    col_idx_ = colIdx;
    n_  = static_cast<int>(rowPtr.size()) - 1;
    nnz_ = static_cast<int>(colIdx.size());
    if (n_ < 0) n_ = 0;

    buildColumnView();
    symbolicFactorization();
    allocateFactors();

    acc_re_.assign(static_cast<size_t>(n_), 0.0);
    acc_im_.assign(static_cast<size_t>(n_), 0.0);
    acc_pos_.assign(static_cast<size_t>(n_), -1);
    acc_gen_.assign(static_cast<size_t>(n_), -1);
    g_work_.assign(static_cast<size_t>(nnz_), 0.0);
    perturbed_ = false;
    min_pivot_ = 0.0;
    last_flops_ = 0;
    message_.clear();
}

void ComplexSparseLU::buildColumnView() {
    // CSR -> CSC (counting sort) plus the mapping back into the CSR values,
    // so a per-frequency assembly reads G/C through value_slot_.
    col_ptr_.assign(static_cast<size_t>(n_) + 1, 0);
    for (int k = 0; k < nnz_; ++k) {
        ++col_ptr_[static_cast<size_t>(col_idx_[k]) + 1];
    }
    for (int j = 0; j < n_; ++j) {
        col_ptr_[static_cast<size_t>(j) + 1] += col_ptr_[static_cast<size_t>(j)];
    }
    row_of_.assign(static_cast<size_t>(nnz_), 0);
    value_slot_.assign(static_cast<size_t>(nnz_), 0);
    std::vector<int> next(col_ptr_.begin(), col_ptr_.end() - 1);
    for (int i = 0; i < n_; ++i) {          // stable: columns end sorted by row
        for (int k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k) {
            const int j = col_idx_[k];
            const int p = next[static_cast<size_t>(j)]++;
            row_of_[static_cast<size_t>(p)] = i;
            value_slot_[static_cast<size_t>(p)] = k;
        }
    }
}

// ===============================================================
// SYMBOLIC FACTORIZATION (once): combined L|U pattern per column
// ===============================================================

void ComplexSparseLU::symbolicFactorization() {
    lcol_ptr_.assign(static_cast<size_t>(n_) + 1, 0);
    lcol_idx_.clear();
    lcol_split_.clear();
    udiag_.clear();
    ll_ = 0;
    if (n_ == 0) return;

    xi_.assign(static_cast<size_t>(n_), 0);
    mark_.assign(static_cast<size_t>(n_), -1);
    stack_.assign(static_cast<size_t>(n_), 0);

    std::vector<int>& Lp = lcol_ptr_;
    std::vector<int> col;      // this column's rows, ascending

    // Column reach (Gilbert-Peierls), no pivoting, natural order: the pattern
    // of column k is the pattern of A(:,k) closed under "for every row j < k
    // it touches, also the pattern of column j". Rows >= k are L entries and
    // need no expansion; j == k is the column itself, so only j < k expands.
    for (int k = 0; k < n_; ++k) {
        int top = 0;
        for (int e = col_ptr_[k]; e < col_ptr_[k + 1]; ++e) {
            xi_[static_cast<size_t>(top++)] = row_of_[static_cast<size_t>(e)];
        }
        for (int p = 0; p < top; ++p) {
            mark_[static_cast<size_t>(xi_[static_cast<size_t>(p)])] = k;
        }
        int sp = 0;
        for (int p = 0; p < top; ++p) {
            const int j = xi_[static_cast<size_t>(p)];
            if (j < k) stack_[static_cast<size_t>(sp++)] = j;
        }
        while (sp > 0) {
            const int j = stack_[static_cast<size_t>(--sp)];
            for (int p = Lp[j]; p < Lp[j + 1]; ++p) {
                const int i = lcol_idx_[static_cast<size_t>(p)];
                if (mark_[static_cast<size_t>(i)] != k) {
                    mark_[static_cast<size_t>(i)] = k;
                    if (i < k) stack_[static_cast<size_t>(sp++)] = i;
                    xi_[static_cast<size_t>(top++)] = i;
                }
            }
        }

        col.assign(xi_.begin(), xi_.begin() + top);
        std::sort(col.begin(), col.end());
        lcol_idx_.insert(lcol_idx_.end(), col.begin(), col.end());
        Lp[static_cast<size_t>(k) + 1] = static_cast<int>(lcol_idx_.size());
    }
    ll_ = static_cast<long long>(lcol_idx_.size());

    // Split at the diagonal: rows < j are U(i,j), rows >= j are L(i,j) with a
    // unit diagonal. The diagonal entry is the FIRST L entry of the column.
    lcol_split_.assign(static_cast<size_t>(n_), 0);
    udiag_.assign(static_cast<size_t>(n_), -1);
    for (int j = 0; j < n_; ++j) {
        const int b = lcol_ptr_[j];
        const int e = lcol_ptr_[j + 1];
        int p = b;
        while (p < e && lcol_idx_[static_cast<size_t>(p)] < j) ++p;
        lcol_split_[static_cast<size_t>(j)] = p;
        for (int q = p; q < e; ++q) {
            if (lcol_idx_[static_cast<size_t>(q)] == j) {
                udiag_[static_cast<size_t>(j)] = q;
                break;
            }
        }
    }
}

void ComplexSparseLU::allocateFactors() {
    const size_t total = static_cast<size_t>(ll_);
    ux_re_.assign(total, 0.0);
    ux_im_.assign(total, 0.0);
    lx_re_.assign(total, 0.0);
    lx_im_.assign(total, 0.0);
}

// ===============================================================
// NUMERIC FACTORIZATION (per frequency)
// ===============================================================

bool ComplexSparseLU::factorize(double omega, const double* G, const double* C) {
    if (n_ == 0) {
        last_flops_ = 0;
        return true;
    }
    perturbed_ = false;
    min_pivot_ = std::numeric_limits<double>::infinity();
    last_flops_ = 0;
    message_.clear();

    std::copy(G, G + nnz_, g_work_.begin());     // mutable copy (perturbation)

    double amax = 0.0;
    for (int k = 0; k < nnz_; ++k) {
        amax = std::max(amax, std::abs(g_work_[static_cast<size_t>(k)]));
        amax = std::max(amax, std::abs(omega * C[k]));
    }
    if (!(amax > 0.0)) {
        message_ = "AC matrix is identically zero";
        return false;
    }
    const double pivot_floor = pivot_tol_ * amax;

    for (int j = 0; j < n_; ++j) {
        const int cb = lcol_ptr_[j];
        const int ce = lcol_ptr_[j + 1];
        const int split = lcol_split_[static_cast<size_t>(j)];
        const int dpos = udiag_[static_cast<size_t>(j)];
        if (dpos < 0) {
            message_ = "AC matrix is structurally singular (row " +
                       std::to_string(j) + ")";
            return false;
        }

        // 1) clear the accumulator on this column's pattern, remember where
        //    each row lives (dense map: every access below is O(1)).
        //    The generation stamp is what makes "is row r part of column j"
        //    an O(1) test WITHOUT resetting the whole map per column: a row
        //    that this column does not touch keeps a stale position, and the
        //    stamp is what tells it apart (reading acc[r] there would inject
        //    the previous column's value).
        for (int p = cb; p < ce; ++p) {
            const int r = lcol_idx_[static_cast<size_t>(p)];
            acc_re_[static_cast<size_t>(r)] = 0.0;
            acc_im_[static_cast<size_t>(r)] = 0.0;
            acc_pos_[static_cast<size_t>(r)] = p;
            acc_gen_[static_cast<size_t>(r)] = j;
        }

        // 2) scatter A(:,j) = G(:,j) + j·omega·C(:,j) through the CSC view.
        for (int e = col_ptr_[j]; e < col_ptr_[j + 1]; ++e) {
            const int r = row_of_[static_cast<size_t>(e)];
            const int slot = value_slot_[static_cast<size_t>(e)];
            if (acc_gen_[static_cast<size_t>(r)] != j) continue;
            acc_re_[static_cast<size_t>(r)] += g_work_[static_cast<size_t>(slot)];
            acc_im_[static_cast<size_t>(r)] += omega * C[slot];
        }

        // 3) for each earlier row i: x[i] is U(i,j); subtract U(i,j)·L(:,i).
        for (int i = 0; i < j; ++i) {
            if (acc_gen_[static_cast<size_t>(i)] != j) continue;  // U(i,j) = 0
            const int ipos = acc_pos_[static_cast<size_t>(i)];
            const double ur = acc_re_[static_cast<size_t>(i)];
            const double ui = acc_im_[static_cast<size_t>(i)];
            ux_re_[static_cast<size_t>(ipos)] = ur;   // U(i,j)
            ux_im_[static_cast<size_t>(ipos)] = ui;
            if (ur == 0.0 && ui == 0.0) continue;
            const int lb = lcol_split_[static_cast<size_t>(i)];
            const int le = lcol_ptr_[i + 1];
            for (int p = lb + 1; p < le; ++p) {  // rows > i (L(i,i) == 1)
                const int r = lcol_idx_[static_cast<size_t>(p)];
                const double lr = lx_re_[static_cast<size_t>(p)];
                const double li = lx_im_[static_cast<size_t>(p)];
                acc_re_[static_cast<size_t>(r)] -= ur * lr - ui * li;
                acc_im_[static_cast<size_t>(r)] -= ur * li + ui * lr;
                last_flops_ += 6;
            }
        }

        // 4) pivot L(j,j) = x[j]; the remaining rows are L(i,j)·L(j,j).
        double pr = acc_re_[static_cast<size_t>(j)];
        double pi = acc_im_[static_cast<size_t>(j)];
        double pmag = std::sqrt(pr * pr + pi * pi);
        if (pmag < pivot_floor) {
            perturbed_ = true;
            const double target = pivot_floor * 10.0;
            int slot = -1;
            for (int e = col_ptr_[j]; e < col_ptr_[j + 1]; ++e) {
                if (row_of_[static_cast<size_t>(e)] == j) {
                    slot = value_slot_[static_cast<size_t>(e)];
                    break;
                }
            }
            if (slot < 0) {
                message_ = "cannot perturb a structurally absent diagonal";
                return false;
            }
            const double wcc = omega * C[slot];
            const double have = std::sqrt(g_work_[static_cast<size_t>(slot)]
                                          * g_work_[static_cast<size_t>(slot)]
                                          + wcc * wcc);
            const double need = std::sqrt(std::max(0.0,
                                   target * target - pi * pi));
            g_work_[static_cast<size_t>(slot)] =
                (need > have) ? need : -(have + target);
            pr = g_work_[static_cast<size_t>(slot)];
            pi = wcc;
            pmag = std::sqrt(pr * pr + pi * pi);
        }
        if (pmag < min_pivot_) min_pivot_ = pmag;
        if (!(pmag > 0.0)) {
            message_ = "AC matrix is singular at column " + std::to_string(j);
            return false;
        }

        // 5) store column j of L (unit diagonal at dpos) and the computed
        //    U(i,j) entries on rows < j. The diagonal slot carries both
        //    factors: U(j,j) = pivot in ux_, L(j,j) = 1 in lx_.
        ux_re_[static_cast<size_t>(dpos)] = pr;
        ux_im_[static_cast<size_t>(dpos)] = pi;
        lx_re_[static_cast<size_t>(dpos)] = 1.0;
        lx_im_[static_cast<size_t>(dpos)] = 0.0;
        const double inv_den = 1.0 / (pr * pr + pi * pi);
        for (int p = split + 1; p < ce; ++p) {
            const int r = lcol_idx_[static_cast<size_t>(p)];
            const double ar = acc_re_[static_cast<size_t>(r)];
            const double ai = acc_im_[static_cast<size_t>(r)];
            lx_re_[static_cast<size_t>(p)] = (ar * pr + ai * pi) * inv_den;
            lx_im_[static_cast<size_t>(p)] = (ai * pr - ar * pi) * inv_den;
        }
        if (std::getenv("EXSPICE_DEBUG_LU")) {
            std::fprintf(stderr, "[LU] col %d: ", j);
            for (int p = cb; p < split; ++p) {
                std::fprintf(stderr, "U(%d,%d)=(%g,%g) ",
                             lcol_idx_[static_cast<size_t>(p)], j,
                             ux_re_[static_cast<size_t>(p)],
                             ux_im_[static_cast<size_t>(p)]);
            }
            std::fprintf(stderr, "| U(%d,%d)=(%g,%g) ", j, j, pr, pi);
            for (int p = split + 1; p < ce; ++p) {
                std::fprintf(stderr, "| L(%d,%d)=(%g,%g) ",
                             lcol_idx_[static_cast<size_t>(p)], j,
                             lx_re_[static_cast<size_t>(p)],
                             lx_im_[static_cast<size_t>(p)]);
            }
            std::fprintf(stderr, "\n");
        }
    }
    return true;
}

// ===============================================================
// TRIANGULAR SOLVES
// ===============================================================

void ComplexSparseLU::solve(std::vector<double>& xr,
                           std::vector<double>& xi) const {
    if (n_ == 0) return;
    // Forward: L y = b (unit diagonal), column order.
    for (int j = 0; j < n_; ++j) {
        const double yr = xr[static_cast<size_t>(j)];
        const double yi = xi[static_cast<size_t>(j)];
        if (yr == 0.0 && yi == 0.0) continue;
        const int lb = lcol_split_[static_cast<size_t>(j)];
        const int le = lcol_ptr_[j + 1];
        for (int p = lb + 1; p < le; ++p) {
            const int r = lcol_idx_[static_cast<size_t>(p)];
            const double lr = lx_re_[static_cast<size_t>(p)];
            const double li = lx_im_[static_cast<size_t>(p)];
            xr[static_cast<size_t>(r)] -= lr * yr - li * yi;
            xi[static_cast<size_t>(r)] -= lr * yi + li * yr;
        }
    }
    // Backward: U x = y, column order. The diagonal division comes FIRST
    // (x[j] must be final before it updates the rows above it).
    for (int j = n_ - 1; j >= 0; --j) {
        const int dpos = udiag_[static_cast<size_t>(j)];
        if (dpos < 0) continue;
        // U(j,j) is the pivot, stored in the diagonal slot of column j.
        const double dr = ux_re_[static_cast<size_t>(dpos)];
        const double di = ux_im_[static_cast<size_t>(dpos)];
        const double den = dr * dr + di * di;
        if (!(den > 0.0)) continue;
        const double nr = xr[static_cast<size_t>(j)];
        const double ni = xi[static_cast<size_t>(j)];
        const double xj = (nr * dr + ni * di) / den;
        const double xij = (ni * dr - nr * di) / den;
        xr[static_cast<size_t>(j)] = xj;
        xi[static_cast<size_t>(j)] = xij;
        if (xj == 0.0 && xij == 0.0) continue;
        const int cb = lcol_ptr_[j];
        const int split = lcol_split_[static_cast<size_t>(j)];
        for (int p = cb; p < split; ++p) {
            const int r = lcol_idx_[static_cast<size_t>(p)];
            const double ur = ux_re_[static_cast<size_t>(p)];
            const double ui = ux_im_[static_cast<size_t>(p)];
            xr[static_cast<size_t>(r)] -= ur * xj - ui * xij;
            xi[static_cast<size_t>(r)] -= ur * xij + ui * xj;
        }
    }
}

} // namespace e_XSpice