/**
 * @author Gabriel Henrique Silva
 * @file MosfetBSIM4Device.cpp
 * @brief BSIM4 level-2 MOSFET — DAE shell (limiters + stamps).
 *
 * Physics lives behind IBsim4Kernel (reference subset by default; Berkeley
 * C via the adapter skeleton when the shim is vendored). All voltages are
 * N-oriented.
 */

#include "MosfetBSIM4Device.h"
#include "../core/DeviceConstants.h"
#include "../external/Bsim4Kernel.h"
#include "../external/Bsim4ReferenceKernel.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace e_XSpice {

bool MosfetBSIM4Device::debug_ = false;

MosfetBSIM4Device::MosfetBSIM4Device(
    RCP<const MosfetBSIM4Model> model,
    double W, double L, double M, double AREA)
    : model_(std::move(model)), kernel_(make_rcp<Bsim4ReferenceKernel>()),
      W_(W), L_(L), M_(M), AREA_(AREA) {
    if (!model_) {
        throw std::invalid_argument("BSIM4 device requires a model card");
    }
    if (W_ <= 0.0 || L_ <= 0.0) {
        throw std::invalid_argument("BSIM4 MOSFET W and L must be positive");
    }
}

double MosfetBSIM4Device::pnjlim(double vnew, double vold, double vt,
                                 double vcrit) {
    if (vnew > vcrit && std::abs(vnew - vold) > 2.0 * vt) {
        if (vold > 0.0) {
            const double arg = 1.0 + (vnew - vold) / vt;
            if (arg > 0.0) {
                vnew = vold + vt * std::log(arg);
            } else {
                vnew = vcrit;
            }
        } else {
            vnew = vt * std::log(vnew / vt);
        }
    }
    return vnew;
}

double MosfetBSIM4Device::clampStep(double vnew, double vold,
                                    double maxStep) {
    const double d = vnew - vold;
    if (d > maxStep) return vold + maxStep;
    if (d < -maxStep) return vold - maxStep;
    return vnew;
}

void MosfetBSIM4Device::limitedVoltages(const DaeSolution& x,
                                         int rd, int rg, int rs, int rb,
                                         bool allow_limiting) const {
    const double vd = (rd >= 0) ? x.data()[rd] : 0.0;
    const double vg = (rg >= 0) ? x.data()[rg] : 0.0;
    const double vs = (rs >= 0) ? x.data()[rs] : 0.0;
    const double vb = (rb >= 0) ? x.data()[rb] : 0.0;

    // N-oriented controlling voltages (p = ±1 applied by caller).
    double vgs = p_ * (vg - vs);
    double vds = p_ * (vd - vs);
    double vbs = p_ * (vb - vs);
    const double vbd = vbs - vds;

    limiting_active_ = false;
    if (allow_limiting && !first_solve_) {
        const double nvgs = clampStep(vgs, vgs_prev_, MaxStep_);
        if (nvgs != vgs) { vgs = nvgs; limiting_active_ = true; }
        const double nvds = clampStep(vds, vds_prev_, MaxStep_);
        if (nvds != vds) { vds = nvds; limiting_active_ = true; }
        const double nvbs = pnjlim(vbs, vbs_prev_, thermal_->vt,
                                   thermal_->vcrit);
        if (nvbs != vbs) { vbs = nvbs; limiting_active_ = true; }
        const double nvbd = pnjlim(vbd, vbd_prev_, thermal_->vt,
                                   thermal_->vcrit);
        if (nvbd != vbd) {
            vbs = nvbd + vds;
            limiting_active_ = true;
        }
    }
    vgs_lim_ = vgs;
    vds_lim_ = vds;
    vbs_lim_ = vbs;

    if (debug_) {
        std::fprintf(stderr,
                     "[MOS %s] lim vgs=%.6g vds=%.6g vbs=%.6g (lim=%d)\n",
                     "MOSFET", vgs, vds, vbs,
                     static_cast<int>(limiting_active_));
    }
}

void MosfetBSIM4Device::evaluate(double vgs, double vds, double vbs) const {
    if (!kernel_) {
        throw std::runtime_error("BSIM4 device has no physics kernel");
    }
    Bsim4KernelInputs in;
    in.card = model_.get();
    in.thermal = thermal_;
    in.W = W_;
    in.L = L_;
    in.M = M_;
    in.AREA = AREA_;
    in.vgs = vgs;
    in.vds = vds;
    in.vbs = vbs;

    Bsim4KernelOutputs out;
    std::string err;
    if (!kernel_->eval(in, out, err)) {
        throw std::runtime_error("BSIM4 kernel '" +
                                 std::string(kernel_->name()) +
                                 "' failed: " + err);
    }

    vgs_cached_ = vgs;
    vds_cached_ = vds;
    vbs_cached_ = vbs;
    vth_cached_ = out.vth;
    saturated_ = out.saturated;
    swapped_ = out.swapped;

    id_ = out.id;
    is_ = out.is;
    ib_ = out.ib;
    gm_ = out.gm;
    gds_ = out.gds;
    gmbs_ = out.gmbs;

    qg_ = out.qg;
    qd_ = out.qd;
    qb_ = out.qb;

    cgs_ = out.cgs;
    cgd_ = out.cgd;
    cgb_ = out.cgb;
    cdg_ = out.cdg;
    cdd_ = out.cdd;
    cdb_ = out.cdb;
    csg_ = out.csg;
    csd_ = out.csd;
    css_ = out.css;
    csb_ = out.csb;
    cbg_ = out.cbg;
    cbd_ = out.cbd;
    cbs_ = out.cbs;
    cbb_ = out.cbb;
}

void MosfetBSIM4Device::loadF(std::vector<double>& F) const {
    // Terminal order: D(0), G(1), S(2), B(3)
    // Terminal currents (positive INTO device)
    if (F.size() < 4) return;
    F[0] += id_;
    F[1] += 0.0;  // gate current = 0
    F[2] += is_;
    F[3] += ib_;
}

void MosfetBSIM4Device::loadQ(std::vector<double>& Q) const {
    if (Q.size() < 4) return;
    // Nodal charges: Qd, Qg, Qs(derived), Qb
    Q[0] += qd_;
    Q[1] += qg_;
    Q[2] += -(qg_ + qd_ + qb_);
    Q[3] += qb_;
}

void MosfetBSIM4Device::loadJacobianF(DaeMatrix& J,
                                      const std::vector<int>& slots) const {
    // slots: 16 entries in row-major order (D,G,S,B x D,G,S,B)
    double ef[16] = {0};

    const double gm = gm_;
    const double gds = gds_;
    const double gmbs = gmbs_;
    const double gibd = 0.0;

    if (!swapped_) {
        ef[0*4+1] = gm;                              // d/dg
        ef[0*4+0] = gds - gibd;                      // d/dd
        ef[0*4+2] = -(gm + gds + gmbs);              // d/ds
        ef[0*4+3] = gmbs + gibd;                     // d/db
        ef[2*4+1] = -gm;                             // s/dg
        ef[2*4+0] = -gds;                            // s/dd
        ef[2*4+2] = (gm + gds + gmbs);               // s/ds
        ef[2*4+3] = -gmbs;                           // s/db
    } else {
        ef[0*4+1] = -gm;                             // d/dg
        ef[0*4+0] = (gm + gds + gmbs) - gibd;        // d/dd
        ef[0*4+2] = -gds;                            // d/ds
        ef[0*4+3] = -gmbs + gibd;                    // d/db
        ef[2*4+1] = gm;                             // s/dg
        ef[2*4+0] = -(gm + gds + gmbs);              // s/dd
        ef[2*4+2] = gds;                            // s/ds
        ef[2*4+3] = gmbs;                           // s/db
    }
    ef[3*4+0] = gibd;                 // b/dd
    ef[3*4+2] = gibd;                 // b/ds
    ef[3*4+3] = -2.0 * gibd;          // b/db

    for (int k = 0; k < 16; ++k) {
        if (slots[k] >= 0) {
            J.accumulateSlot(slots[k], ef[k]);
        }
    }
}

void MosfetBSIM4Device::loadJacobianQ(DaeMatrix& J,
                                      const std::vector<int>& slots) const {
    double eq[16] = {0};

    eq[1*4+1] = cgs_ + cgd_ + cgb_;  // g/dg
    eq[1*4+0] = -cgd_;               // g/dd
    eq[1*4+2] = -cgs_;               // g/ds
    eq[1*4+3] = -cgb_;               // g/db
    eq[0*4+1] = -cgd_;               // d/dg
    eq[0*4+0] = cgd_ + cbd_;         // d/dd
    eq[0*4+3] = -cbd_;               // d/db
    eq[2*4+1] = -cgs_;               // s/dg
    eq[2*4+2] = cgs_ + cbs_;         // s/ds
    eq[2*4+3] = -cbs_;               // s/db
    eq[3*4+1] = -cgb_;               // b/dg
    eq[3*4+0] = -cbd_;               // b/dd
    eq[3*4+2] = -cbs_;               // b/ds
    eq[3*4+3] = cgb_ + cbs_ + cbd_;  // b/db

    for (int k = 0; k < 16; ++k) {
        if (slots[k] >= 0) {
            J.accumulateSlot(slots[k], eq[k]);
        }
    }
}

} // namespace e_XSpice
