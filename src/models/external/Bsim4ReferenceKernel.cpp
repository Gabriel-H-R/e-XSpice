/**
 * @author Gabriel Henrique Silva
 * @file Bsim4ReferenceKernel.cpp
 * @brief Verbatim port of the former MosfetBSIM4Device::evaluate() math.
 */
#include "Bsim4ReferenceKernel.h"
#include "../core_models/Junction.h"
#include "../core/DeviceConstants.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace e_XSpice {

// ε_SiO2 = 3.9·ε_0
static constexpr double EpsOx = 3.9 * 8.854214871e-12;

bool Bsim4ReferenceKernel::eval(const Bsim4KernelInputs& in,
                                Bsim4KernelOutputs& out,
                                std::string& err) const {
    if (!in.card || !in.thermal) {
        err = "reference kernel requires card + thermal snapshots";
        return false;
    }
    if (in.W <= 0.0 || in.L <= 0.0) {
        err = "BSIM4 MOSFET W and L must be positive";
        return false;
    }
    const MosfetBSIM4Model& model = *in.card;
    const MosfetBSIM4Thermal& t = *in.thermal;
    const double mult = in.AREA;
    const double gmin = t.gmin;

    const bool sw = (in.vds < 0.0);
    const double vgs_e = sw ? in.vgs - in.vds : in.vgs;
    const double vds_e = sw ? -in.vds : in.vds;
    const double vbs_e = sw ? in.vbs - in.vds : in.vbs;
    out.swapped = sw;

    const double vto_n = model.VTH0;
    const double sqrtPhi = std::sqrt(std::max(model.PHI, 1e-6));
    const double phib = model.PHI - vbs_e;
    double dvth_dvbs = 0.0;
    double vth = vto_n;
    if (phib > 1e-3) {
        const double sq = std::sqrt(phib);
        vth += model.GAMMA * (sq - sqrtPhi);
        dvth_dvbs = -model.GAMMA / (2.0 * sq);
    } else {
        vth += model.GAMMA * (std::sqrt(1e-3) - sqrtPhi);
    }
    out.vth = vth;

    const double ov = vgs_e - vth;
    const double leff = in.L - 2.0 * model.LD;
    if (leff <= 0.0) {
        err = "BSIM4: Leff = L - 2*LD must be positive";
        return false;
    }
    const double denom_mob = 1.0 + (t.ua + 2.0 * t.ub * ov) * ov / leff;
    const double ueff = t.u0 / denom_mob;

    const double beta = model.KP * (in.W / leff) * in.M;
    const double beta_over_u0 = beta * ueff / t.u0;

    double Id = 0.0, gm = 0.0, gds = 0.0, gmbs = 0.0;
    out.saturated = false;

    if (ov > 0.0) {
        const double vdsat = ov;
        const double lame = 1.0 + model.LAMBDA * vds_e;
        if (vds_e < vdsat) {
            Id = beta_over_u0 * (ov * vds_e - 0.5 * vds_e * vds_e) * lame;
            gm = beta_over_u0 * vds_e * lame;
            gds = beta_over_u0 * (ov - vds_e) * lame +
                  beta_over_u0 * (ov * vds_e - 0.5 * vds_e * vds_e) *
                  model.LAMBDA;
            gmbs = -gm * dvth_dvbs;
        } else {
            Id = 0.5 * beta_over_u0 * ov * ov * lame;
            gm = beta_over_u0 * ov * lame;
            gds = 0.5 * beta_over_u0 * ov * ov * model.LAMBDA;
            gmbs = -gm * dvth_dvbs;
            out.saturated = true;
        }
    }

    const double vbd_p = in.vbs - in.vds;
    const double isBulk = t.isat * mult;
    double arg = std::min(Const::MaxExpArg, in.vbs / t.vt);
    const double evbs = std::exp(arg);
    const double Ibs = isBulk * (evbs - 1.0) + gmin * in.vbs;
    arg = std::min(Const::MaxExpArg, vbd_p / t.vt);
    const double evbd = std::exp(arg);
    const double Ibd = isBulk * (evbd - 1.0) + gmin * vbd_p;

    const double sid = sw ? -1.0 : 1.0;
    out.id = sid * Id + Ibd;
    out.is = -sid * Id + Ibs;
    out.ib = -(Ibs + Ibd);
    out.gm = gm;
    out.gds = gds;
    out.gmbs = gmbs;

    const double cox = EpsOx / model.TOX * in.W * leff * mult;
    const double cgso = model.CGSO * in.W * mult;
    const double cgdo = model.CGDO * in.W * mult;
    const double cgbo = model.CGBO * leff * mult;

    double Cgs_p, Cgd_p, Cgb_p;
    if (!sw) {
        if (out.saturated) {
            Cgs_p = (2.0 / 3.0) * cox + cgso;
            Cgd_p = cgdo;
        } else if (ov > 0.0) {
            Cgs_p = 0.5 * cox + cgso;
            Cgd_p = 0.5 * cox + cgdo;
        } else {
            Cgs_p = cgso;
            Cgd_p = cgdo;
        }
        Cgb_p = cox + cgbo;
    } else {
        if (out.saturated) {
            Cgs_p = cgdo;
            Cgd_p = (2.0 / 3.0) * cox + cgso;
        } else if (ov > 0.0) {
            Cgs_p = cgdo;
            Cgd_p = 0.5 * cox + cgso;
        } else {
            Cgs_p = cgdo;
            Cgd_p = cgso;
        }
        Cgb_p = cox + cgbo;
    }

    const double factor = t.bulkCap.cjo;
    const JunctionGeom geomS{model.CBS * mult * factor, 0.0, model.MJ};
    const JunctionGeom geomD{model.CBD * mult * factor, 0.0, model.MJ};
    const JunctionChargePoint jbs =
        junctionCharge(in.vbs, 0.0, 0.0, geomS, t.bulkCap);
    const JunctionChargePoint jbd =
        junctionCharge(vbd_p, 0.0, 0.0, geomD, t.bulkCap);
    const double Cbs = jbs.c, Cbd = jbd.c;

    const double vgd_p = in.vgs - in.vds;
    const double vgb_p = in.vgs - in.vbs;
    out.qg = Cgs_p * in.vgs + Cgd_p * vgd_p + Cgb_p * vgb_p;
    out.qd = -Cgd_p * vgd_p - jbd.q;
    out.qb = -Cgb_p * vgb_p + jbs.q + jbd.q;

    out.cgs = Cgs_p;
    out.cgd = Cgd_p;
    out.cgb = Cgb_p;
    out.cdg = -Cgd_p;
    out.cdd = Cgd_p;
    out.cdb = -Cbd;
    out.csg = -Cgs_p;
    out.csd = 0.0;
    out.css = Cgs_p;
    out.csb = -Cbs;
    out.cbg = -Cgb_p;
    out.cbd = -Cbd;
    out.cbs = -Cbs;
    out.cbb = Cgb_p;
    return true;
}

} // namespace e_XSpice
