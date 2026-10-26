/**
 * @author Gabriel Henrique Silva
 * @file Mosfet.cpp
 * @brief MOSFET (NMOS and PMOS) — level 1 (Shichman-Hodges) in the shared
 *        DAE formulation.
 */

#include "Mosfet.h"
#include "../core/DaeSolution.h"
#include "../solver/DaeMatrix.h"
#include "../models/core/DeviceConstants.h"
#include "../models/core_models/Junction.h"
#include "../models/core_models/MosfetBSIM4Device.h"
#include "../models/core_models/MosfetBSIM4Model.h"
#include "../models/external/BSim4Adapter.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <stdexcept>

namespace e_XSpice {

using std::exp;
using std::min;
using std::max;

bool Mosfet::debug_ = false;

// Terminal order D(0) G(1) S(2) B(3); pattern index r*4+c.
static inline int midx(int r, int c) { return r * 4 + c; }

// ===============================================================
// Construction
// ===============================================================

Mosfet::Mosfet(const std::string& name,
               Node* drain,
               Node* gate,
               Node* source,
               Node* bulk,
               RCP<const DeviceModel> model,
               const MosfetInstanceParams& inst)
    : Device(name, DeviceType::MOSFET)
    , model_(std::move(model))
    , inst_(inst) {
    if (!model_) {
        throw std::invalid_argument("MOSFET requires a model card");
    }
    if (!drain || !gate || !source || !bulk) {
        throw std::invalid_argument("MOSFET needs 4 terminals (D G S B)");
    }
    if (inst_.W <= 0.0 || inst_.L <= 0.0) {
        throw std::invalid_argument("MOSFET W and L must be positive");
    // Binning range check (Xyce N_DEV_MOSFET_B4.C)
    if (bsim4Model_) {
        if (bsim4Model_->LMIN > 0 && inst_.L < bsim4Model_->LMIN) throw std::invalid_argument("L out of binning range");
        if (bsim4Model_->LMAX > 0 && inst_.L > bsim4Model_->LMAX) throw std::invalid_argument("L out of binning range");
        if (bsim4Model_->WMIN > 0 && inst_.W < bsim4Model_->WMIN) throw std::invalid_argument("W out of binning range");
        if (bsim4Model_->WMAX > 0 && inst_.W > bsim4Model_->WMAX) throw std::invalid_argument("W out of binning range");
    }
    }
    if (inst_.M <= 0.0) {
        throw std::invalid_argument("MOSFET M must be positive");
    }

    // ---- BSIM4 level-2 support ----
    // If the card is LEVEL=2, build a dedicated BSIM4 model card and device
    // evaluator. Two implementations are available:
    //   - Reference C++ (MosfetBSIM4Device): pure C++ subset, always works.
    //   - Berkeley C kernel (BSim4Adapter): wraps the upstream BSIM4 C code
    //     (b4ld.c/b4temp.c) verbatim via the SPICE shim. Selected with
    //     USEBERKELEYKERNEL=1 on the MOSFET instance line.
    // The MosfetModel (level-1) and MosfetBSIM4Model (level-2) are separate
    // classes; we cannot cast between them.
    if (model_->isBSIM4()) {
        bsim4Model_ = rcp_static_cast<const MosfetBSIM4Model>(model_);
        if (!bsim4Model_) {
            throw std::invalid_argument("MOSFET: BSIM4 model card cast failed");
        }
        p_ = bsim4Model_->isNmos() ? 1 : -1;
        useBerkeleyKernel_ = inst_.useBerkeleyKernel;
        p_ = bsim4Model_->isNmos() ? 1 : -1;
        useBerkeleyKernel_ = inst_.useBerkeleyKernel;

        if (useBerkeleyKernel_) {
            // Berkeley C kernel via adapter: wraps BSIM4load/b4temp.c.
            bsim4Adapter_ = new BSim4Adapter(bsim4Model_,
                                             {drain, gate, source, bulk},
                                             inst_.W, inst_.L, inst_.M,
                                             inst_.AD > 0.0 ? inst_.AD : 1.0e-12);
        } else {
            // Reference C++ implementation.
            bsim4Device_ = new MosfetBSIM4Device(bsim4Model_,
                                                 inst_.W, inst_.L, inst_.M,
                                                 inst_.M,  // AREA
                                                 inst_.SA, inst_.SB, inst_.SD,
                                                 inst_.SCA, inst_.SCB, inst_.SCC,
                                                 inst_.SC);
            bsim4Device_->setPolarity(bsim4Model_->isNmos() ? 1 : -1);
            if (inst_.TEMPGiven) {
                ownThermalBSIM4_ = bsim4Model_->computeThermal(inst_.TEMP + inst_.DTEMP);
                bsim4Device_->setThermal(&ownThermalBSIM4_);
            } else {
                bsim4Device_->setThermal(&bsim4Model_->thermal());
            }
        }
    } else {
        p_ = level1Model()->isNmos() ? 1 : -1;
    }

    if (inst_.TEMPGiven) {
        if (inst_.TEMP <= 0.0) {
            throw std::invalid_argument("MOSFET TEMP must be positive (Kelvin)");
        }
        if (!bsim4Model_) {
            ownThermal_ = level1Model()->computeThermal(inst_.TEMP + inst_.DTEMP);
            thermal_ = &ownThermal_;
        }
    } else {
        if (!bsim4Model_) {
            thermal_ = &level1Model()->thermal();
        }
    }

    // Setup-time geometry (never recomputed in the loops).
    // BSIM4 level-2 computes its own geometry internally; level-1 needs
    // these precomputed values. When using the Berkeley C kernel adapter,
    // geometry is handled inside the adapter, so skip level-1 precompute.
    if (!bsim4Device_ && !useBerkeleyKernel_) {
        leff_ = inst_.L - 2.0 * level1Model()->LD;
        if (leff_ <= 0.0) {
            throw std::invalid_argument("MOSFET Leff = L − 2·LD must be positive");
        }
        beta_ = level1Model()->KP * (inst_.W / leff_) * inst_.M;
        cox_  = thermal_->coxPrime * inst_.W * leff_ * inst_.M;
        cgso_ = level1Model()->CGSO * inst_.W * inst_.M;
        cgdo_ = level1Model()->CGDO * inst_.W * inst_.M;
        cgbo_ = level1Model()->CGBO * leff_ * inst_.M;
    }

    // ---- BSIM4 internal-node layout (rgate/rbody/rds) ----
    if (bsim4Model_) {
        int rgatemod = (inst_.RGATEMOD > 0) ? inst_.RGATEMOD : bsim4Model_->RGATEMOD;
        int rdsmod = bsim4Model_->RDSMOD;
        int rbodymod = (inst_.RBODYMOD > 0) ? inst_.RBODYMOD : bsim4Model_->RBODYMOD;
        if (rgatemod > 0) { offG_ = numInternal_++; }
        if (rdsmod > 0) { offD_ = numInternal_++; offS_ = numInternal_++; }
        if (rbodymod > 0) { offB_ = numInternal_++; }
        if (offS_ >= 0) rs_ = (inst_.NRS > 0) ? (bsim4Model_->RDSW * inst_.NRS / (inst_.W * inst_.NF)) : 1e6;
        if (offB_ >= 0) rb_ = (bsim4Model_->RBODYMOD > 0) ? (1.0 / (bsim4Model_->RSH * inst_.NF)) : 1e6;
    }

    nodes_.push_back(drain);   // 0
    nodes_.push_back(gate);    // 1
    nodes_.push_back(source);  // 2
    nodes_.push_back(bulk);    // 3
    drain->addConnectedDevice(this);
    gate->addConnectedDevice(this);
    source->addConnectedDevice(this);
    bulk->addConnectedDevice(this);
}

Mosfet::~Mosfet() {
    delete bsim4Device_;
    delete bsim4Adapter_;
}

// ===============================================================
// CIRCUIT TEMPERATURE (`.TEMP`)
// ===============================================================

void Mosfet::setCircuitTemperature(double kelvin) {
    (void)kelvin;
    // BSIM4 level-2 evaluates its own geometry and thermal quantities
    // inside evaluate() (it reads the card snapshot live), so nothing to
    // refresh. Level-1 precomputed cox at construction time.
    if (bsim4Device_ || inst_.TEMPGiven) return;
    if (bsim4Adapter_) {
        bsim4Adapter_->setCircuitTemperature(kelvin);
        return;
    }
    auto card = level1Model();
    if (!card) return;
    cox_ = thermal_->coxPrime * inst_.W * leff_ * inst_.M;
}

void Mosfet::assignIndices() {
    Device::assignIndices();
    if (bsim4Adapter_) {
        const int r[4] = {rowOf(nodes_[0]), rowOf(nodes_[1]),
                          rowOf(nodes_[2]), rowOf(nodes_[3])};
        bsim4Adapter_->setNodeRows(r[0], r[1], r[2], r[3]);
    }
}

void Mosfet::setJacobianSlots(const std::vector<int>& slots) {
    Device::setJacobianSlots(slots);
    if (bsim4Adapter_) {
        bsim4Adapter_->setJacobianSlots(slots);
    }
}

void Mosfet::setTopology(const std::vector<int>& node_map) {
    Device::setTopology(node_map);
    if (bsim4Adapter_) {
        bsim4Adapter_->setTopology(node_map);
    }
}

RCP<Mosfet> Mosfet::makeNmos(const std::string& name, Node* drain,
                             Node* gate, Node* source, Node* bulk) {
    return make_rcp<Mosfet>(name, drain, gate, source, bulk,
                            MosfetModel::defaultNmos(),
                            MosfetInstanceParams());
}

RCP<Mosfet> Mosfet::makePmos(const std::string& name, Node* drain,
                             Node* gate, Node* source, Node* bulk) {
    return make_rcp<Mosfet>(name, drain, gate, source, bulk,
                            MosfetModel::defaultPmos(),
                            MosfetInstanceParams());
}

// ===============================================================
// Setup
// ===============================================================

int Mosfet::registerStates(StateStore& store) {
    setStateIndex(store.registerStates(3));   // Qg, Qd, Qb
    return 3;
}
void Mosfet::registerJacobianPattern(std::vector<int>& rows,
                                     std::vector<int>& cols) const {
    const int re[4] = {rowOf(nodes_[0]), rowOf(nodes_[1]),
                       rowOf(nodes_[2]), rowOf(nodes_[3])};
    const int ri[4] = {rowDi(), rowGi(), rowSi(), rowBi()};

    // Device Jacobian on intrinsic rows (4x4).
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            rows.push_back(ri[i]);
            cols.push_back(ri[j]);
        }
    }

    // Ohmic resistors (external <-> intrinsic): 4 entries each.
    // Track slot base for each resistor.
    int base = 16;
    if (offD_ >= 0) {
        const_cast<Mosfet*>(this)->patRd_ = base;
        const int pairs[4][2] = {{re[0],re[0]},{re[0],ri[0]},{ri[0],re[0]},{ri[0],ri[0]}};
        for (int k = 0; k < 4; ++k) {
            if (pairs[k][0] >= 0 && pairs[k][1] >= 0) {
                rows.push_back(pairs[k][0]);
                cols.push_back(pairs[k][1]);
                base++;
            }
        }
        const_cast<Mosfet*>(this)->patRd_ = 16; // slot base in jac_slots_
    }
    if (offS_ >= 0) {
        const_cast<Mosfet*>(this)->patRs_ = base;
        const int pairs[4][2] = {{re[2],re[2]},{re[2],ri[2]},{ri[2],re[2]},{ri[2],ri[2]}};
        for (int k = 0; k < 4; ++k) {
            if (pairs[k][0] >= 0 && pairs[k][1] >= 0) {
                rows.push_back(pairs[k][0]);
                cols.push_back(pairs[k][1]);
                base++;
            }
        }
        const_cast<Mosfet*>(this)->patRs_ = (offD_ >= 0) ? 16 + 4 : 16; // approximate
    }
    if (offG_ >= 0) {
        const_cast<Mosfet*>(this)->patRg_ = base;
        const int pairs[4][2] = {{re[1],re[1]},{re[1],ri[1]},{ri[1],re[1]},{ri[1],ri[1]}};
        for (int k = 0; k < 4; ++k) {
            if (pairs[k][0] >= 0 && pairs[k][1] >= 0) {
                rows.push_back(pairs[k][0]);
                cols.push_back(pairs[k][1]);
                base++;
            }
        }
    }
    if (offB_ >= 0) {
        const_cast<Mosfet*>(this)->patRb_ = base;
        const int pairs[4][2] = {{re[3],re[3]},{re[3],ri[3]},{ri[3],re[3]},{ri[3],ri[3]}};
        for (int k = 0; k < 4; ++k) {
            if (pairs[k][0] >= 0 && pairs[k][1] >= 0) {
                rows.push_back(pairs[k][0]);
                cols.push_back(pairs[k][1]);
                base++;
            }
        }
    }
}

void Mosfet::beginSolve(bool initial_solve) {
    // BSIM4 path: the level-1 thermal_ pointer is null for BSIM4 cards
    // (it would segfault below) — the device owns its own history.
    if (bsim4Device_) {
        bsim4Device_->beginSolve(initial_solve);
        if (initial_solve) {
            vbs_prev_ = 0.5;
            vbd_prev_ = 0.0;
            vgs_prev_ = 0.0;
            vds_prev_ = 0.0;
            first_solve_ = true;
        }
        return;
    }
    if (bsim4Adapter_) {
        // Berkeley C kernel adapter handles its own state internally.
        if (initial_solve) {
            vbs_prev_ = 0.5;
            vbd_prev_ = 0.0;
            vgs_prev_ = 0.0;
            vds_prev_ = 0.0;
            first_solve_ = true;
        }
        return;
    }
    if (initial_solve) {
        vbs_prev_ = thermal_->vcrit;
        vbd_prev_ = 0.0;
        vgs_prev_ = 0.0;
        vds_prev_ = 0.0;
        first_solve_ = true;
    }
}

// ===============================================================
// Limiting helpers
// ===============================================================

double Mosfet::pnjlim(double vnew, double vold, double vt, double vcrit) {
    if (vnew > vcrit && std::abs(vnew - vold) > 2.0 * vt) {
        if (vold > 0.0) {
            const double arg = 1.0 + (vnew - vold) / vt;
            if (arg > 0.0) {
                vnew = vold + vt * log(arg);
            } else {
                vnew = vcrit;
            }
        } else {
            vnew = vt * log(vnew / vt);
        }
    }
    return vnew;
}

double Mosfet::clampStep(double vnew, double vold, double maxStep) {
    const double d = vnew - vold;
    if (d > maxStep) return vold + maxStep;
    if (d < -maxStep) return vold - maxStep;
    return vnew;
}

void Mosfet::limitedVoltages(const DaeSolution& x, bool allow_limiting) const {
    const int rd = rowDi();
    const int rg = rowGi();
    const int rs = rowSi();
    const int rb = rowBi();

    const double vd = (rd >= 0) ? x.data()[rd] : 0.0;
    const double vg = (rg >= 0) ? x.data()[rg] : 0.0;
    const double vs = (rs >= 0) ? x.data()[rs] : 0.0;
    const double vb = (rb >= 0) ? x.data()[rb] : 0.0;

    // N-oriented controlling voltages (what the card evaluates).
    double vgs = p_ * (vg - vs);
    double vds = p_ * (vd - vs);
    double vbs = p_ * (vb - vs);
    const double vbd = vbs - vds;

    limiting_active_ = false;
    if (allow_limiting && !first_solve_) {
        // Channel first (bounded step — the SPICE fetlim role): vds must
        // already be sane before the vbd fold below uses it.
        const double nvgs = clampStep(vgs, vgs_prev_, MaxStep_);
        if (nvgs != vgs) { vgs = nvgs; limiting_active_ = true; }
        const double nvds = clampStep(vds, vds_prev_, MaxStep_);
        if (nvds != vds) { vds = nvds; limiting_active_ = true; }
        // Bulk diodes: full exponentials → SPICE3f5 pnjlim.
        const double nvbs = pnjlim(vbs, vbs_prev_, thermal_->vt,
                                   thermal_->vcrit);
        if (nvbs != vbs) { vbs = nvbs; limiting_active_ = true; }
        const double nvbd = pnjlim(vbd, vbd_prev_, thermal_->vt,
                                   thermal_->vcrit);
        if (nvbd != vbd) {
            // Fold the correction back keeping the (already sane) vds.
            vbs = nvbd + vds;
            limiting_active_ = true;
        }
    }
    vgs_lim_ = vgs;
    vds_lim_ = vds;
    vbs_lim_ = vbs;
    }

// ===============================================================
// Evaluation
// ===============================================================

void Mosfet::evaluate(double vgs, double vds, double vbs) const {
    const MosfetThermal& t = *thermal_;
    const double mult = inst_.M;
    const double gmin = t.gmin;

    // ---- source/drain swap for vds < 0 (device is symmetric) ----
    // The effective frame always has vds_e ≥ 0; body effect uses the
    // EFFECTIVE source (vbs_e = vbd_phys when swapped).
    const bool sw = (vds < 0.0);
    swapped_ = sw;
    const double vgs_e = sw ? vgs - vds : vgs;   // vgd when swapped
    const double vds_e = sw ? -vds : vds;
    const double vbs_e = sw ? vbs - vds : vbs;   // vbd when swapped

    // ---- threshold (N-oriented; VTO carries the SPICE sign) ----
    const double vto_n = p_ * level1Model()->VTO;
    const double sqrtPhi = std::sqrt(std::max(level1Model()->PHI, 1e-6));
    const double phib = level1Model()->PHI - vbs_e;   // ≥ PHI in normal operation
    double dvth_dvbs = 0.0;
    double vth = vto_n;
    if (phib > 1e-3) {
        const double sq = std::sqrt(phib);
        vth += level1Model()->GAMMA * (sq - sqrtPhi);
        dvth_dvbs = -level1Model()->GAMMA / (2.0 * sq);
    } else {
        vth += level1Model()->GAMMA * (std::sqrt(1e-3) - sqrtPhi);
    }
    vth_cached_ = vth;

    // ---- channel current + derivatives (effective frame) ----
    const double lam = level1Model()->LAMBDA;
    const double ov = vgs_e - vth;
    double Id = 0.0, gm = 0.0, gds = 0.0, gmbs = 0.0;
    double cgs_e = 0.0, cgd_e = 0.0, cgb_e = 0.0;
    saturated_ = false;
    if (ov > 0.0) {
        const double vdsat = ov;
        if (vds_e < vdsat) {
            // Linear (triode).
            const double lame = 1.0 + lam * vds_e;
            Id = beta_ * (ov * vds_e - 0.5 * vds_e * vds_e) * lame;
            gm = beta_ * vds_e * lame;
            gds = beta_ * (ov - vds_e) * lame +
                  beta_ * (ov * vds_e - 0.5 * vds_e * vds_e) * lam;
            gmbs = -gm * dvth_dvbs;
            cgs_e = 0.5 * cox_;
            cgd_e = 0.5 * cox_;
        } else {
            // Saturation.
            const double lame = 1.0 + lam * vds_e;
            Id = 0.5 * beta_ * ov * ov * lame;
            gm = beta_ * ov * lame;
            gds = 0.5 * beta_ * ov * ov * lam;
            gmbs = -gm * dvth_dvbs;
            cgs_e = (2.0 / 3.0) * cox_;
            saturated_ = true;
        }
    } else {
        // Cutoff: only the gate-bulk capacitance remains.
        cgb_e = cox_;
    }
    gm_ = gm;
    // gds_/gmbs_ below hold the PHYSICAL-frame channel slopes; the stamp
    // uses ef_[] directly, so they are informational only. Store the
    // effective values (documented).
    gds_ = gds;
    gmbs_ = gmbs;

    // ---- bulk diodes (PHYSICAL frame — symmetric, no swap needed) ----
    const double vbd_p = vbs - vds;
    const double isBulk = t.isat * mult;
    double arg = min(Const::MaxExpArg, vbs / t.vt);
    const double evbs = exp(arg);
    const double Ibs = isBulk * (evbs - 1.0) + gmin * vbs;
    const double gibs = isBulk * evbs / t.vt + gmin;
    arg = min(Const::MaxExpArg, vbd_p / t.vt);
    const double evbd = exp(arg);
    const double Ibd = isBulk * (evbd - 1.0) + gmin * vbd_p;
    const double gibd = isBulk * evbd / t.vt + gmin;

    // ---- map channel to the physical frame ----
    // Normal:  id = +Id + Ibd,  is = −Id + Ibs.
    // Swapped: id = −Id + Ibd,  is = +Id + Ibs (Id flows S→D).
    const double sid = sw ? -1.0 : 1.0;
    id_ = sid * Id + Ibd;
    is_ = -sid * Id + Ibs;
    ib_ = -(Ibs + Ibd);

    // Meyer + overlap in the physical frame (swap exchanges Cgs/Cgd).
    const double Cgs_p = (sw ? cgd_e : cgs_e) + cgso_;
    const double Cgd_p = (sw ? cgs_e : cgd_e) + cgdo_;
    const double Cgb_p = cgb_e + cgbo_;

    // Bulk depletion (exact Q/C from the shared junction model).
    const double factor = t.bulkCap.cjo;   // temperature factor of unity cap
    const JunctionGeom geomS{level1Model()->CBS * mult * factor, 0.0, level1Model()->MJ};
    const JunctionGeom geomD{level1Model()->CBD * mult * factor, 0.0, level1Model()->MJ};
    const JunctionChargePoint jbs = junctionCharge(vbs, 0.0, 0.0, geomS,
                                                   t.bulkCap);
    const JunctionChargePoint jbd = junctionCharge(vbd_p, 0.0, 0.0, geomD,
                                                   t.bulkCap);
    const double Cbs = jbs.c, Cbd = jbd.c;

    // ---- nodal charges (charge-conserving by construction) ----
    const double vgd_p = vgs - vds;
    const double vgb_p = vgs - vbs;
    qg_ = Cgs_p * vgs + Cgd_p * vgd_p + Cgb_p * vgb_p;
    qd_ = -Cgd_p * vgd_p - jbd.q;
    // Qs is derived (charge conservation: Qs = −Qg−Qd−Qb, see loadQ).
    qb_ = -Cgb_p * vgb_p + jbs.q + jbd.q;

    // ---- unified physical-frame stamps (row-major D G S B) ----
    // F: rows [id | 0 | is | ib], cols d g s b.
    for (int k = 0; k < 16; ++k) ef_[k] = 0.0;
    for (int k = 0; k < 16; ++k) eq_[k] = 0.0;
    if (!sw) {
        ef_[midx(0, 1)] = gm;                              // d/dg
        ef_[midx(0, 0)] = gds - gibd;                      // d/dd
        ef_[midx(0, 2)] = -(gm + gds + gmbs);              // d/ds
        ef_[midx(0, 3)] = gmbs + gibd;                     // d/db
        ef_[midx(2, 1)] = -gm;                             // s/dg
        ef_[midx(2, 0)] = -gds;                            // s/dd
        ef_[midx(2, 2)] = (gm + gds + gmbs) - gibs;        // s/ds
        ef_[midx(2, 3)] = -gmbs + gibs;                    // s/db
    } else {
        ef_[midx(0, 1)] = -gm;                             // d/dg
        ef_[midx(0, 0)] = (gm + gds + gmbs) - gibd;        // d/dd
        ef_[midx(0, 2)] = -gds;                            // d/ds
        ef_[midx(0, 3)] = -gmbs + gibd;                    // d/db
        ef_[midx(2, 1)] = gm;                              // s/dg
        ef_[midx(2, 0)] = -(gm + gds + gmbs);              // s/dd
        ef_[midx(2, 2)] = gds - gibs;                      // s/ds
        ef_[midx(2, 3)] = gmbs + gibs;                     // s/db
    }
    ef_[midx(3, 0)] = gibd;                 // b/dd
    ef_[midx(3, 2)] = gibs;                 // b/ds
    ef_[midx(3, 3)] = -(gibs + gibd);       // b/db
    // (gate row stays zero: no DC gate current)

    // Q: same layout with the capacitance stamp.
    eq_[midx(1, 1)] = Cgs_p + Cgd_p + Cgb_p;  // g/dg
    eq_[midx(1, 0)] = -Cgd_p;                 // g/dd
    eq_[midx(1, 2)] = -Cgs_p;                 // g/ds
    eq_[midx(1, 3)] = -Cgb_p;                 // g/db
    eq_[midx(0, 1)] = -Cgd_p;                 // d/dg
    eq_[midx(0, 0)] = Cgd_p + Cbd;            // d/dd
    eq_[midx(0, 3)] = -Cbd;                   // d/db
    eq_[midx(2, 1)] = -Cgs_p;                 // s/dg
    eq_[midx(2, 2)] = Cgs_p + Cbs;            // s/ds
    eq_[midx(2, 3)] = -Cbs;                   // s/db
    eq_[midx(3, 1)] = -Cgb_p;                 // b/dg
    eq_[midx(3, 0)] = -Cbd;                   // b/dd
    eq_[midx(3, 2)] = -Cbs;                   // b/ds
    eq_[midx(3, 3)] = Cgb_p + Cbs + Cbd;      // b/db
}

// ===============================================================
// Loaders
// ===============================================================

void Mosfet::pullBSIM4() const {
    id_ = bsim4Device_->getDrainCurrent();
    is_ = bsim4Device_->getSourceCurrent();
    ib_ = bsim4Device_->getBulkCurrent();
    gm_ = bsim4Device_->getGm();
    qg_ = bsim4Device_->getChargeG();
    qd_ = bsim4Device_->getChargeD();
    qb_ = bsim4Device_->getChargeB();
    vth_cached_ = bsim4Device_->getVth();
    saturated_ = bsim4Device_->inSaturation();
    limiting_active_ = bsim4Device_->limitingActive();
}

void Mosfet::loadF(const DaeSolution& x, double time,
                   std::vector<double>& F) const {
    (void)time;
    // BSIM4 level-2 path with reference C++ device.
    if (bsim4Device_) {
        // Intrinsic rows (alias externals when no internal nodes).
        const int ri[4] = {rowDi(), rowGi(), rowSi(), rowBi()};
        const int re[4] = {rowOf(nodes_[0]), rowOf(nodes_[1]), rowOf(nodes_[2]), rowOf(nodes_[3])};
        bsim4Device_->limitedVoltages(x, ri[0], ri[1], ri[2], ri[3], true);
        bsim4Device_->evaluate(bsim4Device_->vgsLim(),
                               bsim4Device_->vdsLim(),
                               bsim4Device_->vbsLim());
        bsim4Device_->loadF(F, ri[0], ri[1], ri[2], ri[3]);
        // Ohmic resistors (external -> intrinsic): physical, no polarity.
        if (false) {
            const double vd_ext = (re[0] >= 0) ? x.data()[re[0]] : 0.0;
            const double vd_int = (ri[0] >= 0) ? x.data()[ri[0]] : 0.0;
            if (re[0] >= 0) F[re[0]] += rd_ * (vd_ext - vd_int);
            if (ri[0] >= 0) F[ri[0]] += rd_ * (vd_int - vd_ext);
        }
        if (offS_ >= 0 && rs_ > 0.0) {
            const double vs_ext = (re[2] >= 0) ? x.data()[re[2]] : 0.0;
            const double vs_int = (ri[2] >= 0) ? x.data()[ri[2]] : 0.0;
            if (re[2] >= 0) F[re[2]] += rs_ * (vs_ext - vs_int);
            if (ri[2] >= 0) F[ri[2]] += rs_ * (vs_int - vs_ext);
        }
        if (offG_ >= 0 && rg_ > 0.0) {
            const double vg_ext = (re[1] >= 0) ? x.data()[re[1]] : 0.0;
            const double vg_int = (ri[1] >= 0) ? x.data()[ri[1]] : 0.0;
            if (re[1] >= 0) F[re[1]] += rg_ * (vg_ext - vg_int);
            if (ri[1] >= 0) F[ri[1]] += rg_ * (vg_int - vg_ext);
        }
        if (offB_ >= 0 && rb_ > 0.0) {
            const double vb_ext = (re[3] >= 0) ? x.data()[re[3]] : 0.0;
            const double vb_int = (ri[3] >= 0) ? x.data()[ri[3]] : 0.0;
            if (re[3] >= 0) F[re[3]] += rb_ * (vb_ext - vb_int);
            if (ri[3] >= 0) F[ri[3]] += rb_ * (vb_int - vb_ext);
        }
        pullBSIM4();
        bsim4Device_->commitLimited();
        first_solve_ = false;
        return;
    }
    // BSIM4 level-2 path with Berkeley C kernel adapter.
    if (bsim4Adapter_) {
        // Set up RHS pointer and call adapter's loadF (which calls
        // BSIM4load internally, stamping ceq directly into F).
        bsim4Adapter_->setRhsPtr(F.data());
        bsim4Adapter_->loadF(x, time, F);
        // Pull outputs from adapter for reporting.
        id_ = bsim4Adapter_->getDrainCurrent();
        is_ = bsim4Adapter_->getSourceCurrent();
        ib_ = bsim4Adapter_->getBulkCurrent();
        gm_ = bsim4Adapter_->getGm();
        vth_cached_ = bsim4Adapter_->getVth();
        saturated_ = false;  // adapter doesn't track saturation in this subset
        limiting_active_ = false;
        first_solve_ = false;
        return;
    }

    limitedVoltages(x, true);
    evaluate(vgs_lim_, vds_lim_, vbs_lim_);
    {
        const int rd = rowOf(nodes_[0]);
        const int rg = rowOf(nodes_[1]);
        const int rs = rowOf(nodes_[2]);
        const int rb = rowOf(nodes_[3]);
        const double vd = (rd >= 0) ? x.data()[rd] : 0.0;
        const double vg = (rg >= 0) ? x.data()[rg] : 0.0;
        const double vs = (rs >= 0) ? x.data()[rs] : 0.0;
        const double vb = (rb >= 0) ? x.data()[rb] : 0.0;
        vgs_cached_ = p_ * (vg - vs);
        vds_cached_ = p_ * (vd - vs);
        vbs_cached_ = p_ * (vb - vs);
    }

    const int r[4] = {rowOf(nodes_[0]), rowOf(nodes_[1]),
                      rowOf(nodes_[2]), rowOf(nodes_[3])};
    const double iv[4] = {id_, 0.0, is_, ib_};
    for (int i = 0; i < 4; ++i) {
        if (r[i] >= 0) F[r[i]] += p_ * iv[i];
    }

    vbs_prev_ = vbs_lim_;
    vbd_prev_ = vbs_lim_ - vds_lim_;
    vgs_prev_ = vgs_lim_;
    vds_prev_ = vds_lim_;
    first_solve_ = false;
}

void Mosfet::applySolutionLimit(std::vector<double>& x) const {
    // BSIM4 path: the device owns the limited voltages and history.
    if (bsim4Device_) {
        bsim4Device_->applySolutionLimit(x, rowDi(), rowGi(),
                                         rowSi(), rowBi(), p_);
        return;
    }
    if (!limiting_active_) return;

    // Rebuild the unknowns from the LIMITED controlling voltages, anchored
    // on the GATE (kept as-is): the gate is the control terminal — usually
    // the stiffest node around — while drain/source/bulk carry the
    // excursion. Moving a voltage-driven gate would fight its source for a
    // whole iteration (observed divergence on a CMOS inverter edge).
    const int rd = rowOf(nodes_[0]);
    const int rg = rowOf(nodes_[1]);
    const int rs = rowOf(nodes_[2]);
    const int rb = rowOf(nodes_[3]);

    const double vg = (rg >= 0) ? x[rg] : 0.0;
    if (rs >= 0) x[rs] = vg - p_ * vgs_lim_;
    const double vs = (rs >= 0) ? x[rs] : 0.0;
    if (rd >= 0) x[rd] = vs + p_ * vds_lim_;
    if (rb >= 0) x[rb] = vs + p_ * vbs_lim_;
}

void Mosfet::loadQ(const DaeSolution& x, double time,
                   std::vector<double>& Q) const {
    (void)time;
    // BSIM4 level-2 path with reference C++ device.
    if (bsim4Device_) {
        const int ri[4] = {rowDi(), rowGi(), rowSi(), rowBi()};
        bsim4Device_->limitedVoltages(x, ri[0], ri[1], ri[2], ri[3], true);
        bsim4Device_->evaluate(bsim4Device_->vgsLim(),
                               bsim4Device_->vdsLim(),
                               bsim4Device_->vbsLim());
        bsim4Device_->loadQ(Q, ri[0], ri[1], ri[2], ri[3]);
        pullBSIM4();
        bsim4Device_->commitLimited();
        first_solve_ = false;
        return;
    }
    // BSIM4 level-2 path with Berkeley C kernel adapter (no charge model).
    if (bsim4Adapter_) {
        first_solve_ = false;
        return;
    }

    limitedVoltages(x, true);
    evaluate(vgs_lim_, vds_lim_, vbs_lim_);

    const int r[4] = {rowOf(nodes_[0]), rowOf(nodes_[1]),
                      rowOf(nodes_[2]), rowOf(nodes_[3])};
    const double qv[4] = {qd_, qg_, -(qg_ + qd_ + qb_), qb_};
    for (int i = 0; i < 4; ++i) {
        if (r[i] >= 0) Q[r[i]] += p_ * qv[i];
    }
}

void Mosfet::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                           double time) const {
    (void)x; (void)time;
    // BSIM4 level-2 path with reference C++ device.
    if (bsim4Device_) {
        // Device Jacobian on intrinsic rows (first 16 slots).
        std::vector<int> devSlots(jac_slots_.begin(), jac_slots_.begin() + 16);
        bsim4Device_->loadJacobianF(J, devSlots);
        // Ohmic resistor conductances.
        int slotBase = 16;
        if (false) {
            // entries: (re,re), (re,ri), (ri,re), (ri,ri)
            if (slotBase + 3 < (int)jac_slots_.size()) {
                J.accumulateSlot(jac_slots_[slotBase+0], rd_);
                J.accumulateSlot(jac_slots_[slotBase+1], -rd_);
                J.accumulateSlot(jac_slots_[slotBase+2], -rd_);
                J.accumulateSlot(jac_slots_[slotBase+3], rd_);
            }
            slotBase += 4;
        }
        if (offS_ >= 0 && rs_ > 0.0) {
            if (slotBase + 3 < (int)jac_slots_.size()) {
                J.accumulateSlot(jac_slots_[slotBase+0], rs_);
                J.accumulateSlot(jac_slots_[slotBase+1], -rs_);
                J.accumulateSlot(jac_slots_[slotBase+2], -rs_);
                J.accumulateSlot(jac_slots_[slotBase+3], rs_);
            }
            slotBase += 4;
        }
        if (offG_ >= 0 && rg_ > 0.0) {
            if (slotBase + 3 < (int)jac_slots_.size()) {
                J.accumulateSlot(jac_slots_[slotBase+0], rg_);
                J.accumulateSlot(jac_slots_[slotBase+1], -rg_);
                J.accumulateSlot(jac_slots_[slotBase+2], -rg_);
                J.accumulateSlot(jac_slots_[slotBase+3], rg_);
            }
            slotBase += 4;
        }
        if (offB_ >= 0 && rb_ > 0.0) {
            if (slotBase + 3 < (int)jac_slots_.size()) {
                J.accumulateSlot(jac_slots_[slotBase+0], rb_);
                J.accumulateSlot(jac_slots_[slotBase+1], -rb_);
                J.accumulateSlot(jac_slots_[slotBase+2], -rb_);
                J.accumulateSlot(jac_slots_[slotBase+3], rb_);
            }
            slotBase += 4;
        }
        return;
    }
    // BSIM4 level-2 path with Berkeley C kernel adapter.
    if (bsim4Adapter_) {
        // Set up matrix values pointer for the SMPaddElt shim, then call
        // adapter's loadJacobianF (which calls BSIM4load internally and
        // emits cached conductances into J via pre-resolved slots).
        bsim4Adapter_->setMatrixValues(J.values().data());
        bsim4Adapter_->loadJacobianF(J, x, time);
        return;
    }
    for (int k = 0; k < 16; ++k) emit(J, k, ef_[k]);
}

void Mosfet::loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                           double time) const {
    (void)x; (void)time;
    // BSIM4 level-2 path with reference C++ device.
    if (bsim4Device_) {
        std::vector<int> devSlots(jac_slots_.begin(), jac_slots_.begin() + 16);
        bsim4Device_->loadJacobianQ(J, devSlots);
        return;
    }
    // BSIM4 level-2 path with Berkeley C kernel adapter (no charge model).
    if (bsim4Adapter_) {
        return;
    }
    for (int k = 0; k < 16; ++k) emit(J, k, eq_[k]);
}

void Mosfet::storeState(const DaeSolution& x, StateStore& store) const {
    // BSIM4 level-2 path with reference C++ device.
    if (bsim4Device_) {
        bsim4Device_->limitedVoltages(x, rowDi(), rowGi(),
                                    rowSi(), rowBi(),
                                    true);
        bsim4Device_->evaluate(bsim4Device_->vgsLim(),
                               bsim4Device_->vdsLim(),
                               bsim4Device_->vbsLim());
        pullBSIM4();
        const double q[3] = {qg_, qd_, qb_};
        store.setDeviceQ(state_index_, q, 3);
        return;
    }
    // BSIM4 level-2 path with Berkeley C kernel adapter (no charge model).
    if (bsim4Adapter_) {
        return;
    }

    limitedVoltages(x, true);
    evaluate(vgs_lim_, vds_lim_, vbs_lim_);
    const double q[3] = {qg_, qd_, qb_};
    store.setDeviceQ(state_index_, q, 3);
}

void Mosfet::loadChargeHistory(const StateStore& store, double q_coeff,
                               double qdot_coeff,
                               std::vector<double>& rhs) const {
    if (q_coeff == 0.0) return;

    // BSIM4 level-2 path with reference C++ device.
    if (bsim4Device_) {
        const int rd = rowDi();
        const int rg = rowGi();
        const int rs = rowSi();
        const int rb = rowBi();

        const double qg = q_coeff * store.deviceQ(state_index_, 0) +
                          qdot_coeff * store.deviceQDot(state_index_, 0);
        const double qd = q_coeff * store.deviceQ(state_index_, 1) +
                          qdot_coeff * store.deviceQDot(state_index_, 1);
        const double qb = q_coeff * store.deviceQ(state_index_, 2) +
                          qdot_coeff * store.deviceQDot(state_index_, 2);
        const double qs = -(qg + qd + qb);

        if (rd >= 0) rhs[rd] += p_ * qd;
        if (rg >= 0) rhs[rg] += p_ * qg;
        if (rs >= 0) rhs[rs] += p_ * qs;
        if (rb >= 0) rhs[rb] += p_ * qb;
        return;
    }
    // BSIM4 level-2 path with Berkeley C kernel adapter (no charge model).
    if (bsim4Adapter_) {
        return;
    }

    const int rd = rowOf(nodes_[0]);
    const int rg = rowOf(nodes_[1]);
    const int rs = rowOf(nodes_[2]);
    const int rb = rowOf(nodes_[3]);

    const double qg = q_coeff * store.deviceQ(state_index_, 0) +
                      qdot_coeff * store.deviceQDot(state_index_, 0);
    const double qd = q_coeff * store.deviceQ(state_index_, 1) +
                      qdot_coeff * store.deviceQDot(state_index_, 1);
    const double qb = q_coeff * store.deviceQ(state_index_, 2) +
                      qdot_coeff * store.deviceQDot(state_index_, 2);
    const double qs = -(qg + qd + qb);

    if (rd >= 0) rhs[rd] += p_ * qd;
    if (rg >= 0) rhs[rg] += p_ * qg;
    if (rs >= 0) rhs[rs] += p_ * qs;
    if (rb >= 0) rhs[rb] += p_ * qb;
}

void Mosfet::refreshOutputs(const DaeSolution& x, double time) const {
    (void)time;
    // BSIM4 level-2 path with reference C++ device (no limiting on converged point).
    if (bsim4Device_) {
        bsim4Device_->limitedVoltages(x, rowDi(), rowGi(),
                                    rowSi(), rowBi(),
                                    false);
        bsim4Device_->evaluate(bsim4Device_->vgsLim(),
                               bsim4Device_->vdsLim(),
                               bsim4Device_->vbsLim());
        pullBSIM4();
        vgs_cached_ = p_ * bsim4Device_->vgsLim();
        vds_cached_ = p_ * bsim4Device_->vdsLim();
        vbs_cached_ = p_ * bsim4Device_->vbsLim();
        return;
    }
    // BSIM4 level-2 path with Berkeley C kernel adapter.
    if (bsim4Adapter_) {
        // Adapter already updated its caches during loadF; just pull outputs.
        id_ = bsim4Adapter_->getDrainCurrent();
        is_ = bsim4Adapter_->getSourceCurrent();
        ib_ = bsim4Adapter_->getBulkCurrent();
        gm_ = bsim4Adapter_->getGm();
        vth_cached_ = bsim4Adapter_->getVth();
        return;
    }

    // Final cache refresh from the converged solution (no limiting).
    const int rd = rowDi();
    const int rg = rowGi();
    const int rs = rowSi();
    const int rb = rowBi();

    const double vd = (rd >= 0) ? x.data()[rd] : 0.0;
    const double vg = (rg >= 0) ? x.data()[rg] : 0.0;
    const double vs = (rs >= 0) ? x.data()[rs] : 0.0;
    const double vb = (rb >= 0) ? x.data()[rb] : 0.0;

    const double vgs = p_ * (vg - vs);
    const double vds = p_ * (vd - vs);
    const double vbs = p_ * (vb - vs);
    vgs_lim_ = vgs;
    vds_lim_ = vds;
    vbs_lim_ = vbs;
    evaluate(vgs, vds, vbs);
    vgs_cached_ = p_ * vgs;
    vds_cached_ = p_ * vds;
    vbs_cached_ = p_ * vbs;
}

bool Mosfet::acceptSolution() const {
    if (bsim4Device_) return bsim4Device_->acceptSolution();
    if (bsim4Adapter_) return bsim4Adapter_->acceptSolution();
    return !limiting_active_;
}

void Mosfet::resetIterationState(const DaeSolution& x) {
    // Re-anchor to the RAW N-oriented voltages at the restored point.
    limiting_active_ = false;
    const int rd = rowOf(nodes_[0]);
    const int rg = rowOf(nodes_[1]);
    const int rs = rowOf(nodes_[2]);
    const int rb = rowOf(nodes_[3]);
    const double vd = (rd >= 0) ? x.data()[rd] : 0.0;
    const double vg = (rg >= 0) ? x.data()[rg] : 0.0;
    const double vs = (rs >= 0) ? x.data()[rs] : 0.0;
    const double vb = (rb >= 0) ? x.data()[rb] : 0.0;
    vbs_prev_ = p_ * (vb - vs);
    const double vds = p_ * (vd - vs);
    vbd_prev_ = vbs_prev_ - vds;
    vgs_prev_ = p_ * (vg - vs);
    vds_prev_ = vds;
    if (bsim4Device_) {
        bsim4Device_->resetHistory(vgs_prev_, vds_prev_, vbs_prev_);
    }
    // Berkeley adapter handles its own history internally.
}

void Mosfet::confirmAccept() {
    if (bsim4Device_) {
        bsim4Device_->confirmAccept();
        limiting_active_ = false;
        return;
    }
    if (bsim4Adapter_) {
        bsim4Adapter_->confirmAccept();
        limiting_active_ = false;
        return;
    }
    limiting_active_ = false;
}

// ===============================================================
// Outputs
// ===============================================================

double Mosfet::getPower() const {
    const double vd = nodes_[0] ? nodes_[0]->getVoltage() : 0.0;
    const double vs = nodes_[2] ? nodes_[2]->getVoltage() : 0.0;
    const double vb = nodes_[3] ? nodes_[3]->getVoltage() : 0.0;
    return vd * (p_ * id_) + vs * (p_ * is_) + vb * (p_ * ib_);
    // (gate draws no DC current; the vg term is identically zero)
}

// ===============================================================
// Noise (`.NOISE`)
// ===============================================================
// Channel thermal noise (Van der Ziel / SPICE3 MOS1): the drain-to-source
// current generator has density (8/3)·k·T·gm. gm_ is the transconductance
// cached by the last evaluation, i.e. exactly the operating-point value that
// the linearization used. gm <= 0 (cutoff) means no channel noise.

void Mosfet::loadNoiseSources(const NoiseContext& ctx,
                              std::vector<NoiseSource>& out) const {
    const double gm = std::fabs(gm_);
    // Thermal noise (drain-source)
    if (gm > 0.0) {
        NoiseSource s;
        s.label = "thermal";
        s.pos = rowOf(nodes_[0]);
        s.neg = rowOf(nodes_[2]);
        s.psd_a2 = (8.0 / 3.0) * ctx.boltzmann * ctx.temperature * gm;
        s.f_exponent = 0.0;
        out.push_back(s);
    }

    // Flicker (1/f) noise for BSIM4
    if (isBSIM4() && bsim4Model_) {
        double Id = std::fabs(getDrainCurrent());
        if (Id > 0.0) {
            NoiseSource f;
            f.label = "flicker";
            f.pos = rowOf(nodes_[0]);
            f.neg = rowOf(nodes_[2]);
            if (bsim4Model_->FNOIMOD == 0) {
                // Simple: KF * Id^AF / (f^EF * Cox * L * W)
                double psd = bsim4Model_->KF * std::pow(Id, bsim4Model_->AF);
                f.psd_a2 = psd;
                f.f_exponent = bsim4Model_->EF;
            } else {
                // BSIM4 flicker noise (NOIA/NOIB/NOIC)
                double Leff = inst_.L - 2.0 * bsim4Model_->LD;
                if (Leff <= 0) Leff = 1e-9;
                double Cox = 3.9 * 8.854e-12 / bsim4Model_->TOX;
                double Nl = Cox * (vgs_cached_ - vth_cached_) / 1.60219e-19;
                if (Nl < 0) Nl = 0;
                double nstar = 1e12;
                double T1 = bsim4Model_->NOIA + bsim4Model_->NOIB / (Nl + nstar)
                          + bsim4Model_->NOIC / ((Nl + nstar) * (Nl + nstar));
                double psd = 1.38065e-23 * ctx.temperature * Id * Id / (Cox * Leff * inst_.W * inst_.M)
                           * T1;
                f.psd_a2 = psd;
                f.f_exponent = 1.0; // 1/f
            }
            out.push_back(f);
        }
    }
}

} // namespace e_XSpice
