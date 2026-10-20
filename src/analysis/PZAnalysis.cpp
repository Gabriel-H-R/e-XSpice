/**
 * @author Gabriel Henrique Silva
 * @file PZAnalysis.cpp
 * @brief Pole-zero analysis.
 */

#include "PZAnalysis.h"

#include "../core/Device.h"
#include "../core/Node.h"
#include "DaeSystem.h"
#include "SmallSignalLinearization.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace e_XSpice {

namespace {
bool iequals(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(a[i])) !=
            std::tolower(static_cast<unsigned char>(b[i]))) return false;
    }
    return true;
}

bool denseSolve(std::vector<std::complex<double>>& A, int n,
                std::vector<std::complex<double>>& b) {
    for (int k = 0; k < n; ++k) {
        int piv = k;
        double best = std::abs(A[static_cast<size_t>(k) * n + static_cast<size_t>(k)]);
        for (int i = k + 1; i < n; ++i) {
            double mag = std::abs(A[static_cast<size_t>(i) * n + static_cast<size_t>(k)]);
            if (mag > best) { best = mag; piv = i; }
        }
        if (best < 1e-30) return false;
        if (piv != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(A[static_cast<size_t>(k) * n + static_cast<size_t>(j)],
                          A[static_cast<size_t>(piv) * n + static_cast<size_t>(j)]);
            }
            std::swap(b[static_cast<size_t>(k)], b[static_cast<size_t>(piv)]);
        }
        std::complex<double> pivval = A[static_cast<size_t>(k) * n + static_cast<size_t>(k)];
        for (int i = k + 1; i < n; ++i) {
            std::complex<double> factor = A[static_cast<size_t>(i) * n + static_cast<size_t>(k)] / pivval;
            for (int j = k; j < n; ++j) {
                A[static_cast<size_t>(i) * n + static_cast<size_t>(j)] -=
                    factor * A[static_cast<size_t>(k) * n + static_cast<size_t>(j)];
            }
            b[static_cast<size_t>(i)] -= factor * b[static_cast<size_t>(k)];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        std::complex<double> sum = b[static_cast<size_t>(i)];
        for (int j = i + 1; j < n; ++j) {
            sum -= A[static_cast<size_t>(i) * n + static_cast<size_t>(j)] *
                   b[static_cast<size_t>(j)];
        }
        b[static_cast<size_t>(i)] = sum / A[static_cast<size_t>(i) * n + static_cast<size_t>(i)];
    }
    return true;
}

std::complex<double> detComplex(const std::vector<std::complex<double>>& A0, int n) {
    std::vector<std::complex<double>> A = A0;
    std::complex<double> det = 1.0;
    for (int k = 0; k < n; ++k) {
        int piv = k;
        double best = std::abs(A[static_cast<size_t>(k) * n + static_cast<size_t>(k)]);
        for (int i = k + 1; i < n; ++i) {
            double mag = std::abs(A[static_cast<size_t>(i) * n + static_cast<size_t>(k)]);
            if (mag > best) { best = mag; piv = i; }
        }
        if (best < 1e-30) return 0.0;
        if (piv != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(A[static_cast<size_t>(k) * n + static_cast<size_t>(j)],
                          A[static_cast<size_t>(piv) * n + static_cast<size_t>(j)]);
            }
            det = -det;
        }
        std::complex<double> pivval = A[static_cast<size_t>(k) * n + static_cast<size_t>(k)];
        det *= pivval;
        for (int i = k + 1; i < n; ++i) {
            std::complex<double> factor = A[static_cast<size_t>(i) * n + static_cast<size_t>(k)] / pivval;
            for (int j = k; j < n; ++j) {
                A[static_cast<size_t>(i) * n + static_cast<size_t>(j)] -=
                    factor * A[static_cast<size_t>(k) * n + static_cast<size_t>(j)];
            }
        }
    }
    return det;
}
} // namespace

PZAnalysis::PZAnalysis(Circuit& circuit) : circuit_(circuit) {}
PZAnalysis::~PZAnalysis() = default;

PZResult PZAnalysis::run() {
    PZResult result;
    result.output_nodes = out1_ + (out2_.empty() || out2_ == "0" ? "" : "," + out2_);
    result.input_nodes = in1_ + (in2_.empty() || in2_ == "0" ? "" : "," + in2_);
    result.current_input = current_input_;

    if (out1_.empty() || in1_.empty()) {
        result.success = false;
        result.message = "missing output or input nodes";
        return result;
    }

    SmallSignalLinearization lin;
    {
        SmallSignalLinearization::Options opt;
        opt.max_iterations = max_iterations_;
        opt.abs_tol = std::max(1e-12, tolerance_ * 1e-3);
        opt.verbose = verbose_;
        opt.use_distributed_layout = false;
        if (!lin.prepare(circuit_, opt)) {
            result.success = false;
            result.message = lin.message();
            return result;
        }
    }
    const int nn = lin.size();
    const auto& node_map = lin.system().nodeMap();

    auto rowOfName = [&](const std::string& name) -> int {
        if (name.empty() || name == "0") return -1;
        for (const auto& n : circuit_.getNodes()) {
            if (n->isGround()) continue;
            if (!iequals(n->getName(), name)) continue;
            const int id = n->getId();
            if (id < 0 || id >= static_cast<int>(node_map.size())) return -1;
            return lin.acRow(node_map[static_cast<size_t>(id)]);
        }
        return -1;
    };

    std::vector<double> c_re(static_cast<size_t>(nn), 0.0);
    std::vector<double> c_im(static_cast<size_t>(nn), 0.0);
    int o1 = rowOfName(out1_), o2 = rowOfName(out2_);
    if (o1 >= 0) c_re[static_cast<size_t>(o1)] += 1.0;
    if (o2 >= 0) c_re[static_cast<size_t>(o2)] -= 1.0;

    std::vector<double> b_re(static_cast<size_t>(nn), 0.0);
    std::vector<double> b_im(static_cast<size_t>(nn), 0.0);
    if (current_input_) {
        int i1 = rowOfName(in1_), i2 = rowOfName(in2_);
        if (i1 >= 0) b_re[static_cast<size_t>(i1)] -= 1.0;
        if (i2 >= 0) b_re[static_cast<size_t>(i2)] += 1.0;
    } else {
        int i1 = rowOfName(in1_), i2 = rowOfName(in2_);
        if (i1 >= 0) b_re[static_cast<size_t>(i1)] -= 1.0;
        if (i2 >= 0) b_re[static_cast<size_t>(i2)] += 1.0;
    }

    const auto& g = lin.gValues();
    const auto& c = lin.cValues();
    const auto& rp = lin.rowPointers();
    const auto& ci = lin.colIndices();

    // ---- DC gain ---------------------------------------------------------
    if (current_input_) {
        std::vector<std::complex<double>> A(
            static_cast<size_t>(nn) * static_cast<size_t>(nn), 0.0);
        std::vector<std::complex<double>> x(static_cast<size_t>(nn));
        for (int i = 0; i < nn; ++i) {
            for (int p = rp[i]; p < rp[i + 1]; ++p) {
                A[static_cast<size_t>(i) * nn + static_cast<size_t>(ci[p])] =
                    std::complex<double>(g[p], 0.0);
            }
            x[static_cast<size_t>(i)] = std::complex<double>(b_re[static_cast<size_t>(i)], 0.0);
        }
        denseSolve(A, nn, x);
        double v_out = 0.0, v_in = 0.0;
        if (o1 >= 0) v_out += x[static_cast<size_t>(o1)].real();
        if (o2 >= 0) v_out -= x[static_cast<size_t>(o2)].real();
        int i1 = rowOfName(in1_), i2 = rowOfName(in2_);
        if (i1 >= 0) v_in += x[static_cast<size_t>(i1)].real();
        if (i2 >= 0) v_in -= x[static_cast<size_t>(i2)].real();
        result.dc_gain = (std::fabs(v_in) > 1e-30)
                             ? std::fabs(v_out / v_in) : 0.0;
    } else {
        // vol mode: augmented system with V(in) = 1.
        int na = nn + 1;
        std::vector<std::complex<double>> Aaug(
            static_cast<size_t>(na) * static_cast<size_t>(na), 0.0);
        std::vector<std::complex<double>> baug(static_cast<size_t>(na), 0.0);
        baug[static_cast<size_t>(na - 1)] = std::complex<double>(1.0, 0.0);
        for (int i = 0; i < nn; ++i) {
            for (int p = rp[i]; p < rp[i + 1]; ++p) {
                Aaug[static_cast<size_t>(i) * na + static_cast<size_t>(ci[p])] =
                    std::complex<double>(g[p], 0.0);
            }
        }
        int in1r = rowOfName(in1_), in2r = rowOfName(in2_);
        if (in1r >= 0) Aaug[static_cast<size_t>(in1r) * na + static_cast<size_t>(na - 1)] += 1.0;
        if (in2r >= 0) Aaug[static_cast<size_t>(in2r) * na + static_cast<size_t>(na - 1)] -= 1.0;
        if (in1r >= 0) Aaug[static_cast<size_t>(na - 1) * na + static_cast<size_t>(in1r)] += 1.0;
        if (in2r >= 0) Aaug[static_cast<size_t>(na - 1) * na + static_cast<size_t>(in2r)] -= 1.0;
        std::vector<std::complex<double>> x = baug;
        bool ok = denseSolve(Aaug, na, x);
        (void)ok;
        double vo = 0.0;
        if (o1 >= 0) vo += x[static_cast<size_t>(o1)].real();
        if (o2 >= 0) vo -= x[static_cast<size_t>(o2)].real();
        result.dc_gain = std::fabs(vo);
    }

    // ---- poles -----------------------------------------------------------
    if (!zeros_only_) {
        if (current_input_) {
            // cur mode: poles from det(G + s C) = 0.
            std::vector<double> grid;
            for (double mag = 1e-9; mag <= 1e9; mag *= 1.5) {
                grid.push_back(-mag);
                grid.push_back(mag);
            }
            std::sort(grid.begin(), grid.end());
            grid.erase(std::unique(grid.begin(), grid.end()), grid.end());

            auto evalDet = [&](double sval) -> double {
                std::vector<std::complex<double>> A(
                    static_cast<size_t>(nn) * static_cast<size_t>(nn), 0.0);
                for (int i = 0; i < nn; ++i) {
                    for (int p = rp[i]; p < rp[i + 1]; ++p) {
                        A[static_cast<size_t>(i) * nn + static_cast<size_t>(ci[p])] +=
                            std::complex<double>(g[p] + sval * c[p], 0.0);
                    }
                }
                return detComplex(A, nn).real();
            };

            std::vector<double> dets(grid.size());
            for (size_t gi = 0; gi < grid.size(); ++gi) dets[gi] = evalDet(grid[gi]);
            for (size_t gi = 0; gi + 1 < grid.size(); ++gi) {
                if (dets[gi] == 0.0) continue;
                if (dets[gi] * dets[gi + 1] >= 0) continue;
                double lo = grid[gi], hi = grid[gi + 1];
                for (int it = 0; it < 100; ++it) {
                    double mid = (lo + hi) / 2.0;
                    double dmid = evalDet(mid);
                    if (dets[gi] * dmid < 0) hi = mid;
                    else lo = mid;
                }
                PZPoint pz;
                pz.real = (lo + hi) / 2.0;
                pz.imag = 0.0;
                pz.is_pole = true;
                result.poles.push_back(pz);
            }
        } else {
            // vol mode: poles from det(A_aug(s)) = 0.
            int in1r = rowOfName(in1_), in2r = rowOfName(in2_);
            std::vector<double> grid;
            for (double mag = 1e-9; mag <= 1e9; mag *= 1.5) {
                grid.push_back(-mag);
                grid.push_back(mag);
            }
            std::sort(grid.begin(), grid.end());
            grid.erase(std::unique(grid.begin(), grid.end()), grid.end());

            auto evalDet = [&](double sval) -> double {
                int na = nn + 1;
                std::vector<std::complex<double>> A2(
                    static_cast<size_t>(na) * static_cast<size_t>(na), 0.0);
                for (int i = 0; i < nn; ++i) {
                    for (int p = rp[i]; p < rp[i + 1]; ++p) {
                        A2[static_cast<size_t>(i) * na + static_cast<size_t>(ci[p])] +=
                            std::complex<double>(g[p] + sval * c[p], 0.0);
                    }
                }
                if (in1r >= 0) A2[static_cast<size_t>(in1r) * na + static_cast<size_t>(na - 1)] += 1.0;
                if (in2r >= 0) A2[static_cast<size_t>(in2r) * na + static_cast<size_t>(na - 1)] -= 1.0;
                if (in1r >= 0) A2[static_cast<size_t>(na - 1) * na + static_cast<size_t>(in1r)] += 1.0;
                if (in2r >= 0) A2[static_cast<size_t>(na - 1) * na + static_cast<size_t>(in2r)] -= 1.0;
                return detComplex(A2, na).real();
            };

            std::vector<double> dets(grid.size());
            for (size_t gi = 0; gi < grid.size(); ++gi) dets[gi] = evalDet(grid[gi]);
            for (size_t gi = 0; gi + 1 < grid.size(); ++gi) {
                if (dets[gi] == 0.0) continue;
                if (dets[gi] * dets[gi + 1] >= 0) continue;
                double lo = grid[gi], hi = grid[gi + 1];
                for (int it = 0; it < 100; ++it) {
                    double mid = (lo + hi) / 2.0;
                    double dmid = evalDet(mid);
                    if (dets[gi] * dmid < 0) hi = mid;
                    else lo = mid;
                }
                PZPoint pz;
                pz.real = (lo + hi) / 2.0;
                pz.imag = 0.0;
                pz.is_pole = true;
                result.poles.push_back(pz);
            }
        }
    }

    // ---- zeros -----------------------------------------------------------
    if (!poles_only_) {
        if (current_input_) {
            // Zeros = roots of det(Gr + s Cr) = 0 (reduced system, output grounded).
            int out_idx = -1;
            for (int i = 0; i < nn; ++i) {
                if (c_re[static_cast<size_t>(i)] != 0.0) { out_idx = i; break; }
            }
            if (out_idx >= 0 && nn > 1) {
                int nr = nn - 1;
                std::vector<double> Gr(static_cast<size_t>(nr) * static_cast<size_t>(nr), 0.0);
                std::vector<double> Cr(static_cast<size_t>(nr) * static_cast<size_t>(nr), 0.0);
                int ri = 0;
                for (int i = 0; i < nn; ++i) {
                    if (i == out_idx) continue;
                    int cj = 0;
                    for (int j = 0; j < nn; ++j) {
                        if (j == out_idx) continue;
                        for (int p = rp[i]; p < rp[i + 1]; ++p) {
                            if (ci[p] == j) {
                                Gr[static_cast<size_t>(ri) * nr + static_cast<size_t>(cj)] = g[p];
                                Cr[static_cast<size_t>(ri) * nr + static_cast<size_t>(cj)] = c[p];
                                break;
                            }
                        }
                        ++cj;
                    }
                    ++ri;
                }
                // Grid search on real axis.
                std::vector<double> grid;
                for (double mag = 1e-9; mag <= 1e9; mag *= 1.5) {
                    grid.push_back(-mag);
                    grid.push_back(mag);
                }
                std::sort(grid.begin(), grid.end());
                grid.erase(std::unique(grid.begin(), grid.end()), grid.end());

                auto evalDet = [&](double sval) -> double {
                    std::vector<std::complex<double>> A(
                        static_cast<size_t>(nr) * static_cast<size_t>(nr), 0.0);
                    for (int i = 0; i < nr; ++i) {
                        for (int j = 0; j < nr; ++j) {
                            A[static_cast<size_t>(i) * nr + static_cast<size_t>(j)] =
                                std::complex<double>(Gr[static_cast<size_t>(i) * nr + static_cast<size_t>(j)], 0.0) +
                                sval * std::complex<double>(Cr[static_cast<size_t>(i) * nr + static_cast<size_t>(j)], 0.0);
                        }
                    }
                    return detComplex(A, nr).real();
                };

                std::vector<double> dets(grid.size());
                for (size_t gi = 0; gi < grid.size(); ++gi) dets[gi] = evalDet(grid[gi]);
                for (size_t gi = 0; gi + 1 < grid.size(); ++gi) {
                    if (dets[gi] == 0.0) continue;
                    if (dets[gi] * dets[gi + 1] >= 0) continue;
                    double lo = grid[gi], hi = grid[gi + 1];
                    for (int it = 0; it < 100; ++it) {
                        double mid = (lo + hi) / 2.0;
                        double dmid = evalDet(mid);
                        if (dets[gi] * dmid < 0) hi = mid;
                        else lo = mid;
                    }
                    PZPoint pz;
                    pz.real = (lo + hi) / 2.0;
                    pz.imag = 0.0;
                    pz.is_pole = false;
                    result.zeros.push_back(pz);
                }
            }
        } else {
            // vol mode: zeros skipped.
        }
    }

    result.success = true;
    std::ostringstream msg;
    msg << result.poles.size() << " pole(s), " << result.zeros.size()
        << " zero(s), DC gain = " << result.dc_gain;
    result.message = msg.str();
    return result;
}

} // namespace e_XSpice