/**
 * @file Diode.cpp
 * @brief Implementação do diodo com modelo de Shockley e PNJLIM
 */

#include "Diode.h"
#include "../solver/SparseMatrix.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace e_XSpice {

using std::string;
using std::exp;
using std::log;
using std::vector;
using std::complex;
using std::invalid_argument;

Diode::Diode(const string& name,
             Node* anode,
             Node* cathode,
             const ModelParams& params)
    : Device(name, DeviceType::DIODE)
    , params_(params)
    , vd_prev_(0.0)
    , first_iter_(true)
{
    if (params_.IS <= 0) {
        throw invalid_argument("Diode IS must be positive");
    }
    if (params_.N <= 0 || params_.N > 5) {
        throw invalid_argument("Diode N must be in (0, 5]");
    }

    nodes_.push_back(anode);
    nodes_.push_back(cathode);
    if (anode)   anode->addConnectedDevice(this);
    if (cathode) cathode->addConnectedDevice(this);
}

// ==================== MODELO FÍSICO ====================

double Diode::thermalVoltage() const {
    constexpr double k = 1.380649e-23;
    constexpr double q = 1.602176634e-19;
    return k * params_.T / q;
}

double Diode::shockleyCurrent(double V) const {
    double nVt = params_.N * thermalVoltage();
    // Limita argumento para evitar overflow de exp()
    // exp(709) ≈ DBL_MAX em IEEE 754 double
    double exponent = V / nVt;
    if (exponent > 80.0) {
        // Aproximação assintótica: I ≈ IS * exp(exponent)
        // Mas usamos 80 como limite prático (exp(80) ≈ 5.5e34)
        exponent = 80.0;
    }
    return params_.IS * (exp(exponent) - 1.0);
}

double Diode::shockleyConductance(double V) const {
    double nVt = params_.N * thermalVoltage();
    double exponent = V / nVt;
    if (exponent > 80.0) exponent = 80.0;
    // G = IS/(nVt) * exp(V/nVt) = (I_D + IS) / (nVt)
    // Forma mais estável numericamente:
    return (params_.IS / nVt) * exp(exponent);
}

// ==================== PNJLIM ====================

double Diode::pnjlim(double v_new, double v_old,
                     double vt, double vcrit) const {
    /*
     * PNJLIM — PN Junction LIMiter
     * Retirado diretamente do SPICE3F5 (pnjlim.c)
     *
     * vcrit = nVt * ln(nVt / (sqrt(2) * IS))
     *
     * Limita a variação de tensão para evitar que o Newton-Raphson
     * oscile entre valores que causam overflow na exponencial.
     */
    if (v_new > vcrit && abs(v_new - v_old) > 2.0 * vt) {
        if (v_old > 0.0) {
            double arg = (v_new - v_old) / vt;
            if (arg > 0.0) {
                v_new = v_old + vt * log(arg + 1.0);
            } else {
                v_new = v_old - vt * log(1.0 - arg);
            }
        } else {
            v_new = vt * log(v_new / vt);
        }
    }
    // Clamp negativo
    constexpr double v_neg_lim = -80.0; // Abaixo disso: I_D ≈ -IS (constante)
    if (v_new < v_neg_lim * vt) {
        v_new = v_neg_lim * vt;
    }
    return v_new;
}

// ==================== STAMP NEWTON-RAPHSON ====================

void Diode::stamp(SparseMatrix& G, vector<double>& b,
                  const vector<int>& node_map) {
    Node* n_anode   = nodes_[0];
    Node* n_cathode = nodes_[1];

    int i = n_anode->isGround()   ? -1 : node_map[n_anode->getId()];
    int j = n_cathode->isGround() ? -1 : node_map[n_cathode->getId()];

    // Tensão atual ânodo-cátodo (lida dos nós)
    double Va = n_anode->getVoltage();
    double Vc = n_cathode->getVoltage();
    double vd = Va - Vc;

    // Parâmetros para PNJLIM
    double nVt   = params_.N * thermalVoltage();
    double vcrit = nVt * log(nVt / (sqrt(2.0) * params_.IS));

    // Aplica PNJLIM se não é a primeira iteração
    double vd_limited = vd;
    if (!first_iter_) {
        vd_limited = pnjlim(vd, vd_prev_, nVt, vcrit);
    }
    first_iter_ = false;
    vd_prev_ = vd_limited;

    // Corrente e condutância no ponto limitado
    double I_D = shockleyCurrent(vd_limited);
    double G_d = shockleyConductance(vd_limited);

    // Corrente equivalente Norton: I_eq = I_D - G_d * V_D
    // (isto implementa o RHS correto para Newton-Raphson)
    double I_eq = I_D - G_d * vd_limited;

    // Stamp da condutância (Jacobiano do diodo)
    if (i >= 0) {
        G.add(i, i, G_d);
        if (j >= 0) G.add(i, j, -G_d);
        // KCL: corrente entra no anodo com -I_eq
        b[i] -= I_eq;
    }
    if (j >= 0) {
        G.add(j, j, G_d);
        if (i >= 0) G.add(j, i, -G_d);
        // KCL: corrente sai do catodo
        b[j] += I_eq;
    }
}

// ==================== STAMP AC ====================

void Diode::stampAC(vector<vector<complex<double>>>& Z,
                    vector<complex<double>>& b_ac,
                    const vector<int>& node_map,
                    double omega) {
    (void)b_ac;
    Node* n_anode   = nodes_[0];
    Node* n_cathode = nodes_[1];

    int i = n_anode->isGround()   ? -1 : node_map[n_anode->getId()];
    int j = n_cathode->isGround() ? -1 : node_map[n_cathode->getId()];

    // Small-signal conductance no ponto DC
    double vd_dc = n_anode->getVoltage() - n_cathode->getVoltage();
    double G_d   = shockleyConductance(vd_dc);

    // Capacitância de junção (para representação correta em AC)
    double C_j = getJunctionCapacitance();

    // Y_diode = G_d + jω*C_j
    complex<double> y(G_d, omega * C_j);

    if (i >= 0) { Z[i][i] += y; if (j >= 0) Z[i][j] -= y; }
    if (j >= 0) { Z[j][j] += y; if (i >= 0) Z[j][i] -= y; }
}

// ==================== GRANDEZAS CALCULADAS ====================

double Diode::getCurrent() const {
    double vd = nodes_[0]->getVoltage() - nodes_[1]->getVoltage();
    return shockleyCurrent(vd);
}

double Diode::getPower() const {
    return getVoltage() * getCurrent();
}

double Diode::getVoltage() const {
    return nodes_[0]->getVoltage() - nodes_[1]->getVoltage();
}

double Diode::getSmallSignalConductance() const {
    double vd = nodes_[0]->getVoltage() - nodes_[1]->getVoltage();
    return shockleyConductance(vd);
}

double Diode::getJunctionCapacitance() const {
    if (params_.CJO <= 0.0) return 0.0;

    double vd = nodes_[0]->getVoltage() - nodes_[1]->getVoltage();

    // Modelo de capacitância de junção abrupta:
    // C_j(V) = CJ0 / sqrt(1 - V/VJ)   para V < FC*VJ
    // Aproximação linear para V próximo de VJ (evita singularidade)
    constexpr double FC = 0.5; // Coeficiente de depletação
    double VJ = params_.VJ;
    double CJ0 = params_.CJO;
    double M   = params_.M;

    if (vd <= FC * VJ) {
        return CJ0 / pow(1.0 - vd / VJ, M);
    } else {
        // Linearização além de FC*VJ
        double f1 = pow(1.0 - FC, -(1.0 + M));
        double f2 = 1.0 - FC * (1.0 + M) + M * vd / VJ;
        return CJ0 * f1 * f2;
    }
}

} // namespace e_XSpice