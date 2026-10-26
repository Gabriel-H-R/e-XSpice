/**
 * @author Gabriel Henrique Silva
 * @file MosfetBSIM4Model.cpp
 * @brief BSIM4 4.8.3 model card — level 2 MOSFET (NMOS + PMOS).
 */

#include "MosfetBSIM4Model.h"
#include "../core/DeviceConstants.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace e_XSpice {

MosfetBSIM4Model::MosfetBSIM4Model(std::string name, bool nmos)
    : DeviceModel(std::move(name)), nmos_(nmos) {
    // Card parameter table (canonical first, aliases after — canonical wins
    // when a card spells both). Adding a parameter = one line here.
    declareParams({
        makeIntParam("LEVEL", &MosfetBSIM4Model::LEVEL, "2=BSIM4"),
        makeParam("VTH0", &MosfetBSIM4Model::VTH0, "zero-bias threshold voltage"),
        makeParam("KP", &MosfetBSIM4Model::KP, "transconductance parameter"),
        makeParam("LAMBDA", &MosfetBSIM4Model::LAMBDA,
                  "channel-length modulation"),
        makeParam("RDSW", &MosfetBSIM4Model::RDSW,
                  "source/drain series resistance"),
        makeParam("PRWG", &MosfetBSIM4Model::PRWG,
                  "Rds gate-bias coefficient"),
        makeParam("PRWB", &MosfetBSIM4Model::PRWB,
                  "Rds body-bias coefficient"),
        makeParam("PRT", &MosfetBSIM4Model::PRT,
                  "Rds temperature coefficient"),
        makeParam("RBDB", &MosfetBSIM4Model::RBDB,
                  "body-drain resistance"),
        makeParam("RBSB", &MosfetBSIM4Model::RBSB,
                  "body-source resistance"),
        makeParam("RBPB", &MosfetBSIM4Model::RBPB,
                  "body-bulk resistance"),
        makeParam("RBPS", &MosfetBSIM4Model::RBPS,
                  "body-source (prime) resistance"),
        makeParam("RBPD", &MosfetBSIM4Model::RBPD,
                  "body-drain (prime) resistance"),
        makeParam("XRCRG1", &MosfetBSIM4Model::XRCRG1,
                  "gate resistance coefficient 1"),
        makeParam("XRCRG2", &MosfetBSIM4Model::XRCRG2,
                  "gate resistance coefficient 2"),
        makeParam("PHI", &MosfetBSIM4Model::PHI, "surface potential"),
        makeParam("GAMMA", &MosfetBSIM4Model::GAMMA,
                  "bulk threshold parameter"),
        makeParam("K2", &MosfetBSIM4Model::K2,
                  "body-effect coefficient 2"),
        makeParam("CBD", &MosfetBSIM4Model::CBD,
                  "bulk-drain capacitance (total)"),
        makeParam("CBS", &MosfetBSIM4Model::CBS,
                  "bulk-source capacitance (total)"),
        makeParam("IS", &MosfetBSIM4Model::IS, "bulk saturation current"),
        makeParam("PB", &MosfetBSIM4Model::PB, "bulk junction potential"),
        makeParam("MJ", &MosfetBSIM4Model::MJ, "bulk junction grading"),
        makeParam("FC", &MosfetBSIM4Model::FC,
                  "bulk forward-bias coefficient"),
        makeParam("CGSO", &MosfetBSIM4Model::CGSO,
                  "gate-source overlap / width"),
        makeParam("CGDO", &MosfetBSIM4Model::CGDO,
                  "gate-drain overlap / width"),
        makeParam("CGBO", &MosfetBSIM4Model::CGBO,
                  "gate-bulk overlap / length"),
        makeParam("TOX", &MosfetBSIM4Model::TOX, "oxide thickness (m)"),
        makeParam("LD", &MosfetBSIM4Model::LD, "lateral diffusion (m)"),
        makeParam("EG", &MosfetBSIM4Model::EG, "bandgap voltage"),
        makeParam("XTI", &MosfetBSIM4Model::XTI,
                  "saturation-current temperature exponent"),
        makeParam("GMIN", &MosfetBSIM4Model::GMIN, "parallel conductance"),
        // BSIM4-specific
        makeParam("U0", &MosfetBSIM4Model::U0, "zero-bias mobility"),
        makeParam("VSAT", &MosfetBSIM4Model::VSAT, "saturation velocity"),
        makeParam("UA", &MosfetBSIM4Model::UA,
                  "1st-order mobility degradation"),
        makeParam("UB", &MosfetBSIM4Model::UB,
                  "2nd-order mobility degradation"),
        makeParam("UC", &MosfetBSIM4Model::UC,
                  "saturation velocity coefficient"),
        makeParam("UD", &MosfetBSIM4Model::UD, "saturation effect coeff."),
        makeParam("A0", &MosfetBSIM4Model::A0, "bulk charge coeff. 1"),
        makeParam("B0", &MosfetBSIM4Model::B0, "bulk charge coeff. 2"),
        makeParam("B1", &MosfetBSIM4Model::B1, "bulk charge coeff. 3"),
        makeParam("KETA", &MosfetBSIM4Model::KETA,
                  "back-bias threshold coeff."),
        makeParam("KT1", &MosfetBSIM4Model::KT1,
                  "threshold temp coeff. 1"),
        makeParam("KT2", &MosfetBSIM4Model::KT2,
                  "threshold temp coeff. 2"),
        makeParam("ALPHA0", &MosfetBSIM4Model::ALPHA0,
                  "substrate current coefficient"),
        makeParam("AGIDL", &MosfetBSIM4Model::AGIDL,
                  "GIDL coefficient"),
        makeParam("BGIDL", &MosfetBSIM4Model::BGIDL,
                  "GIDL exponential coefficient"),
        makeParam("CGIDL", &MosfetBSIM4Model::CGIDL,
                  "GIDL body-effect coefficient"),
        makeParam("EGIDL", &MosfetBSIM4Model::EGIDL,
                  "GIDL activation energy"),
        makeParam("AIGC", &MosfetBSIM4Model::AIGC,
                  "gate current coeff. accumulation"),
        makeParam("BIGC", &MosfetBSIM4Model::BIGC,
                  "gate current coeff. accumulation"),
        makeParam("CIGC", &MosfetBSIM4Model::CIGC,
                  "gate current coeff. accumulation"),
        makeParam("AIGS", &MosfetBSIM4Model::AIGS,
                  "gate current coeff. source"),
        makeParam("AIGD", &MosfetBSIM4Model::AIGD,
                  "gate current coeff. drain"),
        makeParam("AIGBACC", &MosfetBSIM4Model::AIGBACC,
                  "bulk current coeff. accumulation"),
        makeParam("AIGBINV", &MosfetBSIM4Model::AIGBINV,
                  "bulk current coeff. inversion"),
        makeParam("NIGC", &MosfetBSIM4Model::NIGC,
                  "gate current exponent"),
        makeIntParam("TRNQSMOD", &MosfetBSIM4Model::TRNQSMOD,
                     "transient NQS model"),
        makeIntParam("ACNQSMOD", &MosfetBSIM4Model::ACNQSMOD,
                     "AC NQS model"),
        makeParam("TAUN", &MosfetBSIM4Model::TAUN,
                  "NQS time constant"),
        makeParam("LMIN", &MosfetBSIM4Model::LMIN,
                  "minimum channel length for binning"),
        makeParam("LMAX", &MosfetBSIM4Model::LMAX,
                  "maximum channel length for binning"),
        makeParam("WMIN", &MosfetBSIM4Model::WMIN,
                  "minimum channel width for binning"),
        makeParam("WMAX", &MosfetBSIM4Model::WMAX,
                  "maximum channel width for binning"),
        makeParam("DMCG", &MosfetBSIM4Model::DMCG,
                  "distance gate to contact"),
        makeParam("SA", &MosfetBSIM4Model::SA,
                  "distance to STI (source)"),
        makeParam("SB", &MosfetBSIM4Model::SB,
                  "distance to STI (drain)"),
        makeParam("SD", &MosfetBSIM4Model::SD,
                  "distance to STI (end)"),
        makeParam("SCA", &MosfetBSIM4Model::SCA,
                  "stress (source)"),
        makeParam("SCB", &MosfetBSIM4Model::SCB,
                  "stress (drain)"),
        makeParam("SCC", &MosfetBSIM4Model::SCC,
                  "stress (channel)"),
        makeParam("SC", &MosfetBSIM4Model::SC,
                  "stress (combined)"),
        makeParam("KVTH0WE", &MosfetBSIM4Model::KVTH0WE,
                  "well proximity Vth coeff."),
        makeParam("K2WE", &MosfetBSIM4Model::K2WE,
                  "well proximity K2 coeff."),
        makeParam("KU0WE", &MosfetBSIM4Model::KU0WE,
                  "well proximity mobility coeff."),
        makeParam("DMCI", &MosfetBSIM4Model::DMCI,
                  "distance gate to contact (inner)"),
        makeParam("DMDG", &MosfetBSIM4Model::DMDG,
                  "distance gate to contact (drain)"),
        makeParam("KF", &MosfetBSIM4Model::KF,
                  "flicker noise coefficient"),
        makeParam("NFACTOR", &MosfetBSIM4Model::NFACTOR,
                  "subthreshold ideality factor"),
        makeParam("VOFF", &MosfetBSIM4Model::VOFF,
                  "threshold voltage offset"),
        makeParam("AF", &MosfetBSIM4Model::AF,
                  "flicker noise exponent"),
        makeParam("EF", &MosfetBSIM4Model::EF,
                  "flicker noise frequency exponent"),
        makeParam("RNOIA", &MosfetBSIM4Model::RNOIA,
                  "thermal noise coeff. drain"),
        makeParam("RNOIB", &MosfetBSIM4Model::RNOIB,
                  "thermal noise coeff. source"),
        makeParam("RNOIC", &MosfetBSIM4Model::RNOIC,
                  "thermal noise coeff. gate"),
        makeParam("AGISL", &MosfetBSIM4Model::AGISL,
                  "GISL coefficient"),
        makeParam("BGISL", &MosfetBSIM4Model::BGISL,
                  "GISL exponential coefficient"),
        makeParam("CGISL", &MosfetBSIM4Model::CGISL,
                  "GISL body-effect coefficient"),
        makeParam("EGISL", &MosfetBSIM4Model::EGISL,
                  "GISL activation energy"),
        makeParam("ALPHA1", &MosfetBSIM4Model::ALPHA1,
                  "substrate current exponential coefficient"),
        makeParam("UTE", &MosfetBSIM4Model::UTE,
                  "mobility temperature exponent"),
        makeParam("DVT0", &MosfetBSIM4Model::DVT0,
                  "short-channel threshold coefficient 1"),
        makeParam("DVT1", &MosfetBSIM4Model::DVT1,
                  "short-channel threshold coefficient 2"),
        makeParam("DVT2", &MosfetBSIM4Model::DVT2,
                  "short-channel threshold coefficient 3"),
        makeParam("DVT2W", &MosfetBSIM4Model::DVT2W,
                  "short-channel threshold coeff. width"),
        makeParam("DVTP0", &MosfetBSIM4Model::DVTP0,
                  "pocket implant threshold coeff."),
        makeParam("W0", &MosfetBSIM4Model::W0, "narrow width offset"),
        makeParam("NFS", &MosfetBSIM4Model::NFS,
                  "subthreshold slope coefficient (optional)"),
        // Xyce BSIM4 4.8.2 addPar (N_DEV_MOSFET_B4.C), non-binned
        makeParam("EOT", &MosfetBSIM4Model::EOT, "Xyce B4.C:339"),
        makeParam("VDDEOT", &MosfetBSIM4Model::VDDEOT, "Xyce B4.C:345"),
        makeParam("TEMPEOT", &MosfetBSIM4Model::TEMPEOT, "Xyce B4.C:350"),
        makeParam("LEFFEOT", &MosfetBSIM4Model::LEFFEOT, "Xyce B4.C:356"),
        makeParam("WEFFEOT", &MosfetBSIM4Model::WEFFEOT, "Xyce B4.C:362"),
        makeParam("ADOS", &MosfetBSIM4Model::ADOS, "Xyce B4.C:368"),
        makeParam("BDOS", &MosfetBSIM4Model::BDOS, "Xyce B4.C:373"),
        makeParam("TOXE", &MosfetBSIM4Model::TOXE, "Xyce B4.C:378"),
        makeParam("TOXP", &MosfetBSIM4Model::TOXP, "Xyce B4.C:384"),
        makeParam("TOXM", &MosfetBSIM4Model::TOXM, "Xyce B4.C:390"),
        makeParam("TOXREF", &MosfetBSIM4Model::TOXREF, "Xyce B4.C:396"),
        makeParam("DTOX", &MosfetBSIM4Model::DTOX, "Xyce B4.C:401"),
        makeParam("EPSROX", &MosfetBSIM4Model::EPSROX, "Xyce B4.C:407"),
        makeParam("CDSC", &MosfetBSIM4Model::CDSC, "Xyce B4.C:413"),
        makeParam("CDSCB", &MosfetBSIM4Model::CDSCB, "Xyce B4.C:418"),
        makeParam("CDSCD", &MosfetBSIM4Model::CDSCD, "Xyce B4.C:423"),
        makeParam("CIT", &MosfetBSIM4Model::CIT, "Xyce B4.C:428"),
        makeParam("XJ", &MosfetBSIM4Model::XJ, "Xyce B4.C:438"),
        makeParam("AT", &MosfetBSIM4Model::AT, "Xyce B4.C:448"),
        makeParam("AGS", &MosfetBSIM4Model::AGS, "Xyce B4.C:458"),
        makeParam("A1", &MosfetBSIM4Model::A1, "Xyce B4.C:463"),
        makeParam("A2", &MosfetBSIM4Model::A2, "Xyce B4.C:468"),
        makeParam("PHIG", &MosfetBSIM4Model::PHIG, "Xyce B4.C:478"),
        makeParam("EPSRGATE", &MosfetBSIM4Model::EPSRGATE, "Xyce B4.C:484"),
        makeParam("EASUB", &MosfetBSIM4Model::EASUB, "Xyce B4.C:489"),
        makeParam("EPSRSUB", &MosfetBSIM4Model::EPSRSUB, "Xyce B4.C:494"),
        makeParam("NI0SUB", &MosfetBSIM4Model::NI0SUB, "Xyce B4.C:499"),
        makeParam("BG0SUB", &MosfetBSIM4Model::BG0SUB, "Xyce B4.C:504"),
        makeParam("TBGASUB", &MosfetBSIM4Model::TBGASUB, "Xyce B4.C:509"),
        makeParam("TBGBSUB", &MosfetBSIM4Model::TBGBSUB, "Xyce B4.C:515"),
        makeParam("NSUB", &MosfetBSIM4Model::NSUB, "Xyce B4.C:520"),
        makeParam("NDEP", &MosfetBSIM4Model::NDEP, "Xyce B4.C:526"),
        makeParam("NSD", &MosfetBSIM4Model::NSD, "Xyce B4.C:532"),
        makeParam("CKAPPAS", &MosfetBSIM4Model::CKAPPAS, "Xyce B4.C:1416"),
        makeParam("CKAPPAD", &MosfetBSIM4Model::CKAPPAD, "Xyce B4.C:1421"),
        makeParam("CF", &MosfetBSIM4Model::CF, "Xyce B4.C:1479"),
        makeParam("CLC", &MosfetBSIM4Model::CLC, "Xyce B4.C:1432"),
        makeParam("CLE", &MosfetBSIM4Model::CLE, "Xyce B4.C:1437"),
        makeParam("ACDE", &MosfetBSIM4Model::ACDE, "Xyce B4.C:1141"),
        makeParam("VFBCV", &MosfetBSIM4Model::VFBCV, "Xyce B4.C:1100"),
        makeParam("VOFFCV", &MosfetBSIM4Model::VOFFCV, "Xyce B4.C:1156"),
        makeParam("VOFFCVL", &MosfetBSIM4Model::VOFFCVL, "Xyce B4.C:806"),
        makeParam("MINVCV", &MosfetBSIM4Model::MINVCV, "Xyce B4.C:796"),
        makeParam("DLC", &MosfetBSIM4Model::DLC, "Xyce B4.C:1447"),
        makeParam("DWC", &MosfetBSIM4Model::DWC, "Xyce B4.C:1442"),
        makeParam("CGDL", &MosfetBSIM4Model::CGDL, "Xyce B4.C:1411"),
        makeParam("CGSL", &MosfetBSIM4Model::CGSL, "Xyce B4.C:1406"),
        makeParam("PHIN", &MosfetBSIM4Model::PHIN, "Xyce B4.C:537"),
        makeParam("NGATE", &MosfetBSIM4Model::NGATE, "Xyce B4.C:542"),
        makeParam("GAMMA1", &MosfetBSIM4Model::GAMMA1, "Xyce B4.C:547"),
        makeParam("GAMMA2", &MosfetBSIM4Model::GAMMA2, "Xyce B4.C:553"),
        makeParam("VBX", &MosfetBSIM4Model::VBX, "Xyce B4.C:559"),
        makeParam("VBM", &MosfetBSIM4Model::VBM, "Xyce B4.C:565"),
        makeParam("XT", &MosfetBSIM4Model::XT, "Xyce B4.C:570"),
        makeParam("K1", &MosfetBSIM4Model::K1, "Xyce B4.C:576"),
        makeParam("KT1L", &MosfetBSIM4Model::KT1L, "Xyce B4.C:587"),
        makeParam("K3", &MosfetBSIM4Model::K3, "Xyce B4.C:603"),
        makeParam("K3B", &MosfetBSIM4Model::K3B, "Xyce B4.C:608"),
        makeParam("DVTP1", &MosfetBSIM4Model::DVTP1, "Xyce B4.C:623"),
        makeParam("DVTP2", &MosfetBSIM4Model::DVTP2, "Xyce B4.C:628"),
        makeParam("DVTP3", &MosfetBSIM4Model::DVTP3, "Xyce B4.C:634"),
        makeParam("DVTP4", &MosfetBSIM4Model::DVTP4, "Xyce B4.C:640"),
        makeParam("DVTP5", &MosfetBSIM4Model::DVTP5, "Xyce B4.C:646"),
        makeParam("LPE0", &MosfetBSIM4Model::LPE0, "Xyce B4.C:652"),
        makeParam("LPEB", &MosfetBSIM4Model::LPEB, "Xyce B4.C:657"),
        makeParam("DVT0W", &MosfetBSIM4Model::DVT0W, "Xyce B4.C:678"),
        makeParam("DVT1W", &MosfetBSIM4Model::DVT1W, "Xyce B4.C:683"),
        makeParam("DROUT", &MosfetBSIM4Model::DROUT, "Xyce B4.C:693"),
        makeParam("DSUB", &MosfetBSIM4Model::DSUB, "Xyce B4.C:698"),
        makeParam("UA1", &MosfetBSIM4Model::UA1, "Xyce B4.C:714"),
        makeParam("UB1", &MosfetBSIM4Model::UB1, "Xyce B4.C:724"),
        makeParam("UC1", &MosfetBSIM4Model::UC1, "Xyce B4.C:734"),
        makeParam("UD1", &MosfetBSIM4Model::UD1, "Xyce B4.C:744"),
        makeParam("UP", &MosfetBSIM4Model::UP, "Xyce B4.C:749"),
        makeParam("LP", &MosfetBSIM4Model::LP, "Xyce B4.C:754"),
        makeParam("EU", &MosfetBSIM4Model::EU, "Xyce B4.C:764"),
        makeParam("UCS", &MosfetBSIM4Model::UCS, "Xyce B4.C:769"),
        makeParam("UCSTE", &MosfetBSIM4Model::UCSTE, "Xyce B4.C:780"),
        makeParam("MINV", &MosfetBSIM4Model::MINV, "Xyce B4.C:791"),
        makeParam("VOFFL", &MosfetBSIM4Model::VOFFL, "Xyce B4.C:801"),
        makeParam("XPART", &MosfetBSIM4Model::XPART, "Xyce B4.C:834"),
        makeParam("DELTA", &MosfetBSIM4Model::DELTA, "Xyce B4.C:839"),
        makeParam("RSH", &MosfetBSIM4Model::RSH, "Xyce B4.C:844"),
        makeParam("RDSWMIN", &MosfetBSIM4Model::RDSWMIN, "Xyce B4.C:854"),
        makeParam("RSW", &MosfetBSIM4Model::RSW, "Xyce B4.C:859"),
        makeParam("RDW", &MosfetBSIM4Model::RDW, "Xyce B4.C:864"),
        makeParam("RDWMIN", &MosfetBSIM4Model::RDWMIN, "Xyce B4.C:869"),
        makeParam("RSWMIN", &MosfetBSIM4Model::RSWMIN, "Xyce B4.C:874"),
        makeParam("ETA0", &MosfetBSIM4Model::ETA0, "Xyce B4.C:894"),
        makeParam("ETAB", &MosfetBSIM4Model::ETAB, "Xyce B4.C:899"),
        makeParam("PCLM", &MosfetBSIM4Model::PCLM, "Xyce B4.C:904"),
        makeParam("PDIBLC1", &MosfetBSIM4Model::PDIBLC1, "Xyce B4.C:910"),
        makeParam("PDIBLC2", &MosfetBSIM4Model::PDIBLC2, "Xyce B4.C:915"),
        makeParam("PDIBLCB", &MosfetBSIM4Model::PDIBLCB, "Xyce B4.C:920"),
        makeParam("FPROUT", &MosfetBSIM4Model::FPROUT, "Xyce B4.C:925"),
        makeParam("PDITS", &MosfetBSIM4Model::PDITS, "Xyce B4.C:930"),
        makeParam("PDITSL", &MosfetBSIM4Model::PDITSL, "Xyce B4.C:935"),
        makeParam("PDITSD", &MosfetBSIM4Model::PDITSD, "Xyce B4.C:940"),
        makeParam("PSCBE1", &MosfetBSIM4Model::PSCBE1, "Xyce B4.C:945"),
        makeParam("PSCBE2", &MosfetBSIM4Model::PSCBE2, "Xyce B4.C:950"),
        makeParam("PVAG", &MosfetBSIM4Model::PVAG, "Xyce B4.C:955"),
        makeParam("JSS", &MosfetBSIM4Model::JSS, "Xyce B4.C:960"),
        makeParam("JSWS", &MosfetBSIM4Model::JSWS, "Xyce B4.C:965"),
        makeParam("JSWGS", &MosfetBSIM4Model::JSWGS, "Xyce B4.C:970"),
        makeParam("PBS", &MosfetBSIM4Model::PBS, "Xyce B4.C:975"),
        makeParam("NJS", &MosfetBSIM4Model::NJS, "Xyce B4.C:980"),
        makeParam("XTIS", &MosfetBSIM4Model::XTIS, "Xyce B4.C:985"),
        makeParam("MJS", &MosfetBSIM4Model::MJS, "Xyce B4.C:990"),
        makeParam("PBSWS", &MosfetBSIM4Model::PBSWS, "Xyce B4.C:995"),
        makeParam("MJSWS", &MosfetBSIM4Model::MJSWS, "Xyce B4.C:1000"),
        makeParam("PBSWGS", &MosfetBSIM4Model::PBSWGS, "Xyce B4.C:1005"),
        makeParam("MJSWGS", &MosfetBSIM4Model::MJSWGS, "Xyce B4.C:1010"),
        makeParam("CJS", &MosfetBSIM4Model::CJS, "Xyce B4.C:1015"),
        makeParam("CJSWS", &MosfetBSIM4Model::CJSWS, "Xyce B4.C:1020"),
        makeParam("CJSWGS", &MosfetBSIM4Model::CJSWGS, "Xyce B4.C:1025"),
        makeParam("JSD", &MosfetBSIM4Model::JSD, "Xyce B4.C:1030"),
        makeParam("JSWD", &MosfetBSIM4Model::JSWD, "Xyce B4.C:1035"),
        makeParam("JSWGD", &MosfetBSIM4Model::JSWGD, "Xyce B4.C:1040"),
        makeParam("PBD", &MosfetBSIM4Model::PBD, "Xyce B4.C:1045"),
        makeParam("NJD", &MosfetBSIM4Model::NJD, "Xyce B4.C:1050"),
        makeParam("XTID", &MosfetBSIM4Model::XTID, "Xyce B4.C:1055"),
        makeParam("MJD", &MosfetBSIM4Model::MJD, "Xyce B4.C:1060"),
        makeParam("PBSWD", &MosfetBSIM4Model::PBSWD, "Xyce B4.C:1065"),
        makeParam("MJSWD", &MosfetBSIM4Model::MJSWD, "Xyce B4.C:1070"),
        makeParam("PBSWGD", &MosfetBSIM4Model::PBSWGD, "Xyce B4.C:1075"),
        makeParam("MJSWGD", &MosfetBSIM4Model::MJSWGD, "Xyce B4.C:1080"),
        makeParam("CJD", &MosfetBSIM4Model::CJD, "Xyce B4.C:1085"),
        makeParam("CJSWD", &MosfetBSIM4Model::CJSWD, "Xyce B4.C:1090"),
        makeParam("CJSWGD", &MosfetBSIM4Model::CJSWGD, "Xyce B4.C:1095"),
        makeParam("TPB", &MosfetBSIM4Model::TPB, "Xyce B4.C:1111"),
        makeParam("TCJ", &MosfetBSIM4Model::TCJ, "Xyce B4.C:1116"),
        makeParam("TPBSW", &MosfetBSIM4Model::TPBSW, "Xyce B4.C:1121"),
        makeParam("TCJSW", &MosfetBSIM4Model::TCJSW, "Xyce B4.C:1126"),
        makeParam("TPBSWG", &MosfetBSIM4Model::TPBSWG, "Xyce B4.C:1131"),
        makeParam("TCJSWG", &MosfetBSIM4Model::TCJSWG, "Xyce B4.C:1136"),
        makeParam("MOIN", &MosfetBSIM4Model::MOIN, "Xyce B4.C:1146"),
        makeParam("NOFF", &MosfetBSIM4Model::NOFF, "Xyce B4.C:1151"),
        makeParam("XGW", &MosfetBSIM4Model::XGW, "Xyce B4.C:1181"),
        makeParam("XGL", &MosfetBSIM4Model::XGL, "Xyce B4.C:1186"),
        makeParam("RSHG", &MosfetBSIM4Model::RSHG, "Xyce B4.C:1191"),
        makeParam("NGCON", &MosfetBSIM4Model::NGCON, "Xyce B4.C:1196"),
        makeParam("VTL", &MosfetBSIM4Model::VTL, "Xyce B4.C:1217"),
        makeParam("LC", &MosfetBSIM4Model::LC, "Xyce B4.C:1223"),
        makeParam("XN", &MosfetBSIM4Model::XN, "Xyce B4.C:1228"),
        makeParam("VFBSDOFF", &MosfetBSIM4Model::VFBSDOFF, "Xyce B4.C:1233"),
        makeParam("TVFBSDOFF", &MosfetBSIM4Model::TVFBSDOFF, "Xyce B4.C:1238"),
        makeParam("TVOFF", &MosfetBSIM4Model::TVOFF, "Xyce B4.C:1243"),
        makeParam("TNFACTOR", &MosfetBSIM4Model::TNFACTOR, "Xyce B4.C:1248"),
        makeParam("TETA0", &MosfetBSIM4Model::TETA0, "Xyce B4.C:1254"),
        makeParam("LINT", &MosfetBSIM4Model::LINT, "Xyce B4.C:1271"),
        makeParam("LL", &MosfetBSIM4Model::LL, "Xyce B4.C:1276"),
        makeParam("LLN", &MosfetBSIM4Model::LLN, "Xyce B4.C:1286"),
        makeParam("LW", &MosfetBSIM4Model::LW, "Xyce B4.C:1291"),
        makeParam("LWN", &MosfetBSIM4Model::LWN, "Xyce B4.C:1301"),
        makeParam("LWL", &MosfetBSIM4Model::LWL, "Xyce B4.C:1306"),
        makeParam("WR", &MosfetBSIM4Model::WR, "Xyce B4.C:1326"),
        makeParam("WINT", &MosfetBSIM4Model::WINT, "Xyce B4.C:1331"),
        makeParam("DWG", &MosfetBSIM4Model::DWG, "Xyce B4.C:1336"),
        makeParam("DWB", &MosfetBSIM4Model::DWB, "Xyce B4.C:1341"),
        makeParam("WL", &MosfetBSIM4Model::WL, "Xyce B4.C:1346"),
        makeParam("WLN", &MosfetBSIM4Model::WLN, "Xyce B4.C:1356"),
        makeParam("WW", &MosfetBSIM4Model::WW, "Xyce B4.C:1361"),
        makeParam("WWN", &MosfetBSIM4Model::WWN, "Xyce B4.C:1371"),
        makeParam("WWL", &MosfetBSIM4Model::WWL, "Xyce B4.C:1376"),
        makeParam("XW", &MosfetBSIM4Model::XW, "Xyce B4.C:1453"),
        makeParam("XL", &MosfetBSIM4Model::XL, "Xyce B4.C:1458"),
        makeParam("BETA0", &MosfetBSIM4Model::BETA0, "Xyce B4.C:1488"),
        makeParam("RGIDL", &MosfetBSIM4Model::RGIDL, "Xyce B4.C:1508"),
        makeParam("KGIDL", &MosfetBSIM4Model::KGIDL, "Xyce B4.C:1514"),
        makeParam("FGIDL", &MosfetBSIM4Model::FGIDL, "Xyce B4.C:1520"),
        makeParam("RGISL", &MosfetBSIM4Model::RGISL, "Xyce B4.C:1546"),
        makeParam("KGISL", &MosfetBSIM4Model::KGISL, "Xyce B4.C:1552"),
        makeParam("FGISL", &MosfetBSIM4Model::FGISL, "Xyce B4.C:1558"),
        makeParam("AIGSD", &MosfetBSIM4Model::AIGSD, "Xyce B4.C:1585"),
        makeParam("BIGSD", &MosfetBSIM4Model::BIGSD, "Xyce B4.C:1590"),
        makeParam("CIGSD", &MosfetBSIM4Model::CIGSD, "Xyce B4.C:1595"),
        makeParam("BIGS", &MosfetBSIM4Model::BIGS, "Xyce B4.C:1605"),
        makeParam("CIGS", &MosfetBSIM4Model::CIGS, "Xyce B4.C:1610"),
        makeParam("BIGD", &MosfetBSIM4Model::BIGD, "Xyce B4.C:1620"),
        makeParam("CIGD", &MosfetBSIM4Model::CIGD, "Xyce B4.C:1625"),
        makeParam("BIGBACC", &MosfetBSIM4Model::BIGBACC, "Xyce B4.C:1635"),
        makeParam("CIGBACC", &MosfetBSIM4Model::CIGBACC, "Xyce B4.C:1640"),
        makeParam("BIGBINV", &MosfetBSIM4Model::BIGBINV, "Xyce B4.C:1650"),
        makeParam("CIGBINV", &MosfetBSIM4Model::CIGBINV, "Xyce B4.C:1655"),
        makeParam("NIGBINV", &MosfetBSIM4Model::NIGBINV, "Xyce B4.C:1665"),
        makeParam("NIGBACC", &MosfetBSIM4Model::NIGBACC, "Xyce B4.C:1670"),
        makeParam("NTOX", &MosfetBSIM4Model::NTOX, "Xyce B4.C:1675"),
        makeParam("EIGBINV", &MosfetBSIM4Model::EIGBINV, "Xyce B4.C:1680"),
        makeParam("PIGCD", &MosfetBSIM4Model::PIGCD, "Xyce B4.C:1685"),
        makeParam("POXEDGE", &MosfetBSIM4Model::POXEDGE, "Xyce B4.C:1691"),
        makeParam("IJTHDFWD", &MosfetBSIM4Model::IJTHDFWD, "Xyce B4.C:1697"),
        makeParam("IJTHSFWD", &MosfetBSIM4Model::IJTHSFWD, "Xyce B4.C:1702"),
        makeParam("IJTHDREV", &MosfetBSIM4Model::IJTHDREV, "Xyce B4.C:1707"),
        makeParam("IJTHSREV", &MosfetBSIM4Model::IJTHSREV, "Xyce B4.C:1712"),
        makeParam("XJBVD", &MosfetBSIM4Model::XJBVD, "Xyce B4.C:1717"),
        makeParam("XJBVS", &MosfetBSIM4Model::XJBVS, "Xyce B4.C:1722"),
        makeParam("BVD", &MosfetBSIM4Model::BVD, "Xyce B4.C:1727"),
        makeParam("BVS", &MosfetBSIM4Model::BVS, "Xyce B4.C:1732"),
        makeParam("JTSS", &MosfetBSIM4Model::JTSS, "Xyce B4.C:1737"),
        makeParam("JTSD", &MosfetBSIM4Model::JTSD, "Xyce B4.C:1742"),
        makeParam("JTSSWS", &MosfetBSIM4Model::JTSSWS, "Xyce B4.C:1747"),
        makeParam("JTSSWD", &MosfetBSIM4Model::JTSSWD, "Xyce B4.C:1752"),
        makeParam("JTSSWGS", &MosfetBSIM4Model::JTSSWGS, "Xyce B4.C:1757"),
        makeParam("JTSSWGD", &MosfetBSIM4Model::JTSSWGD, "Xyce B4.C:1762"),
        makeParam("JTWEFF", &MosfetBSIM4Model::JTWEFF, "Xyce B4.C:1767"),
        makeParam("NJTS", &MosfetBSIM4Model::NJTS, "Xyce B4.C:1773"),
        makeParam("NJTSSW", &MosfetBSIM4Model::NJTSSW, "Xyce B4.C:1778"),
        makeParam("NJTSSWG", &MosfetBSIM4Model::NJTSSWG, "Xyce B4.C:1783"),
        makeParam("NJTSD", &MosfetBSIM4Model::NJTSD, "Xyce B4.C:1788"),
        makeParam("NJTSSWD", &MosfetBSIM4Model::NJTSSWD, "Xyce B4.C:1793"),
        makeParam("NJTSSWGD", &MosfetBSIM4Model::NJTSSWGD, "Xyce B4.C:1798"),
        makeParam("XTSS", &MosfetBSIM4Model::XTSS, "Xyce B4.C:1803"),
        makeParam("XTSD", &MosfetBSIM4Model::XTSD, "Xyce B4.C:1808"),
        makeParam("XTSSWS", &MosfetBSIM4Model::XTSSWS, "Xyce B4.C:1813"),
        makeParam("XTSSWD", &MosfetBSIM4Model::XTSSWD, "Xyce B4.C:1818"),
        makeParam("XTSSWGS", &MosfetBSIM4Model::XTSSWGS, "Xyce B4.C:1823"),
        makeParam("XTSSWGD", &MosfetBSIM4Model::XTSSWGD, "Xyce B4.C:1828"),
        makeParam("TNJTS", &MosfetBSIM4Model::TNJTS, "Xyce B4.C:1833"),
        makeParam("TNJTSSW", &MosfetBSIM4Model::TNJTSSW, "Xyce B4.C:1838"),
        makeParam("TNJTSSWG", &MosfetBSIM4Model::TNJTSSWG, "Xyce B4.C:1843"),
        makeParam("TNJTSD", &MosfetBSIM4Model::TNJTSD, "Xyce B4.C:1848"),
        makeParam("TNJTSSWD", &MosfetBSIM4Model::TNJTSSWD, "Xyce B4.C:1853"),
        makeParam("TNJTSSWGD", &MosfetBSIM4Model::TNJTSSWGD, "Xyce B4.C:1858"),
        makeParam("VTSS", &MosfetBSIM4Model::VTSS, "Xyce B4.C:1863"),
        makeParam("VTSD", &MosfetBSIM4Model::VTSD, "Xyce B4.C:1868"),
        makeParam("VTSSWS", &MosfetBSIM4Model::VTSSWS, "Xyce B4.C:1873"),
        makeParam("VTSSWD", &MosfetBSIM4Model::VTSSWD, "Xyce B4.C:1878"),
        makeParam("VTSSWGS", &MosfetBSIM4Model::VTSSWGS, "Xyce B4.C:1883"),
        makeParam("VTSSWGD", &MosfetBSIM4Model::VTSSWGD, "Xyce B4.C:1888"),
        makeParam("GBMIN", &MosfetBSIM4Model::GBMIN, "Xyce B4.C:1893"),
        makeParam("RBPS0", &MosfetBSIM4Model::RBPS0, "Xyce B4.C:1923"),
        makeParam("RBPSL", &MosfetBSIM4Model::RBPSL, "Xyce B4.C:1929"),
        makeParam("RBPSW", &MosfetBSIM4Model::RBPSW, "Xyce B4.C:1934"),
        makeParam("RBPSNF", &MosfetBSIM4Model::RBPSNF, "Xyce B4.C:1939"),
        makeParam("RBPD0", &MosfetBSIM4Model::RBPD0, "Xyce B4.C:1944"),
        makeParam("RBPDL", &MosfetBSIM4Model::RBPDL, "Xyce B4.C:1950"),
        makeParam("RBPDW", &MosfetBSIM4Model::RBPDW, "Xyce B4.C:1955"),
        makeParam("RBPDNF", &MosfetBSIM4Model::RBPDNF, "Xyce B4.C:1960"),
        makeParam("RBPBX0", &MosfetBSIM4Model::RBPBX0, "Xyce B4.C:1965"),
        makeParam("RBPBXL", &MosfetBSIM4Model::RBPBXL, "Xyce B4.C:1970"),
        makeParam("RBPBXW", &MosfetBSIM4Model::RBPBXW, "Xyce B4.C:1975"),
        makeParam("RBPBXNF", &MosfetBSIM4Model::RBPBXNF, "Xyce B4.C:1980"),
        makeParam("RBPBY0", &MosfetBSIM4Model::RBPBY0, "Xyce B4.C:1985"),
        makeParam("RBPBYL", &MosfetBSIM4Model::RBPBYL, "Xyce B4.C:1990"),
        makeParam("RBPBYW", &MosfetBSIM4Model::RBPBYW, "Xyce B4.C:1995"),
        makeParam("RBPBYNF", &MosfetBSIM4Model::RBPBYNF, "Xyce B4.C:2000"),
        makeParam("RBSBX0", &MosfetBSIM4Model::RBSBX0, "Xyce B4.C:2005"),
        makeParam("RBSBY0", &MosfetBSIM4Model::RBSBY0, "Xyce B4.C:2011"),
        makeParam("RBDBX0", &MosfetBSIM4Model::RBDBX0, "Xyce B4.C:2017"),
        makeParam("RBDBY0", &MosfetBSIM4Model::RBDBY0, "Xyce B4.C:2023"),
        makeParam("RBSDBXL", &MosfetBSIM4Model::RBSDBXL, "Xyce B4.C:2029"),
        makeParam("RBSDBXW", &MosfetBSIM4Model::RBSDBXW, "Xyce B4.C:2034"),
        makeParam("RBSDBXNF", &MosfetBSIM4Model::RBSDBXNF, "Xyce B4.C:2039"),
        makeParam("RBSDBYL", &MosfetBSIM4Model::RBSDBYL, "Xyce B4.C:2044"),
        makeParam("RBSDBYW", &MosfetBSIM4Model::RBSDBYW, "Xyce B4.C:2049"),
        makeParam("RBSDBYNF", &MosfetBSIM4Model::RBSDBYNF, "Xyce B4.C:2054"),
        makeParam("SAREF", &MosfetBSIM4Model::SAREF, "Xyce B4.C:4446"),
        makeParam("SBREF", &MosfetBSIM4Model::SBREF, "Xyce B4.C:4451"),
        makeParam("WLOD", &MosfetBSIM4Model::WLOD, "Xyce B4.C:4456"),
        makeParam("KU0", &MosfetBSIM4Model::KU0, "Xyce B4.C:4461"),
        makeParam("KVSAT", &MosfetBSIM4Model::KVSAT, "Xyce B4.C:4466"),
        makeParam("KVTH0", &MosfetBSIM4Model::KVTH0, "Xyce B4.C:4471"),
        makeParam("TKU0", &MosfetBSIM4Model::TKU0, "Xyce B4.C:4476"),
        makeParam("LLODKU0", &MosfetBSIM4Model::LLODKU0, "Xyce B4.C:4481"),
        makeParam("WLODKU0", &MosfetBSIM4Model::WLODKU0, "Xyce B4.C:4486"),
        makeParam("LLODVTH", &MosfetBSIM4Model::LLODVTH, "Xyce B4.C:4491"),
        makeParam("WLODVTH", &MosfetBSIM4Model::WLODVTH, "Xyce B4.C:4496"),
        makeParam("STK2", &MosfetBSIM4Model::STK2, "Xyce B4.C:4531"),
        makeParam("LODK2", &MosfetBSIM4Model::LODK2, "Xyce B4.C:4536"),
        makeParam("STETA0", &MosfetBSIM4Model::STETA0, "Xyce B4.C:4541"),
        makeParam("LODETA0", &MosfetBSIM4Model::LODETA0, "Xyce B4.C:4546"),
        makeParam("WEB", &MosfetBSIM4Model::WEB, "Xyce B4.C:4552"),
        makeParam("WEC", &MosfetBSIM4Model::WEC, "Xyce B4.C:4557"),
        makeParam("SCREF", &MosfetBSIM4Model::SCREF, "Xyce B4.C:4577"),
        makeParam("GIDLCLAMP", &MosfetBSIM4Model::GIDLCLAMP, "Xyce B4.C:4674"),
        makeParam("IDOVVDS", &MosfetBSIM4Model::IDOVVDS, "Xyce B4.C:4683"),
        makeIntParam("CVCHARGEMOD", &MosfetBSIM4Model::CVCHARGEMOD, "Xyce B4.C:4720"),
        makeIntParam("CAPMOD", &MosfetBSIM4Model::CAPMOD, "Xyce B4.C:4725"),
        makeIntParam("DIOMOD", &MosfetBSIM4Model::DIOMOD, "Xyce B4.C:4730"),
        makeIntParam("RDSMOD", &MosfetBSIM4Model::RDSMOD, "Xyce B4.C:4735"),
        makeIntParam("MOBMOD", &MosfetBSIM4Model::MOBMOD, "Xyce B4.C:4752"),
        makeIntParam("PERMOD", &MosfetBSIM4Model::PERMOD, "Xyce B4.C:4767"),
        makeIntParam("GEOMOD", &MosfetBSIM4Model::GEOMOD, "Xyce B4.C:4772"),
        makeIntParam("RGEOMOD", &MosfetBSIM4Model::RGEOMOD, "Xyce B4.C:4777"),
        makeIntParam("FNOIMOD", &MosfetBSIM4Model::FNOIMOD, "Xyce B4.C:4782"),
        makeIntParam("TNOIMOD", &MosfetBSIM4Model::TNOIMOD, "Xyce B4.C:4787"),
        makeIntParam("MTRLMOD", &MosfetBSIM4Model::MTRLMOD, "Xyce B4.C:4792"),
        makeIntParam("MTRLCOMPATMOD", &MosfetBSIM4Model::MTRLCOMPATMOD, "Xyce B4.C:4797"),
        makeIntParam("TEMPMOD", &MosfetBSIM4Model::TEMPMOD, "Xyce B4.C:4813"),
        makeIntParam("GIDLMOD", &MosfetBSIM4Model::GIDLMOD, "Xyce B4.C:4818"),
        makeIntParam("PARAMCHK", &MosfetBSIM4Model::PARAMCHK, "Xyce B4.C:4824"),
        makeIntParam("BINUNIT", &MosfetBSIM4Model::BINUNIT, "Xyce B4.C:4829"),
        makeParam("DMCGT", &MosfetBSIM4Model::DMCGT, "DMCG temperature coeff."),
        makeIntParam("RBODYMOD", &MosfetBSIM4Model::RBODYMOD, "body resistance model"),
        makeIntParam("RGATEMOD", &MosfetBSIM4Model::RGATEMOD, "gate resistance model"),
        makeIntParam("WPEMOD", &MosfetBSIM4Model::WPEMOD, "well proximity model"),
        makeIntParam("IGCMOD", &MosfetBSIM4Model::IGCMod, "gate current model"),
        makeIntParam("IGBMOD", &MosfetBSIM4Model::IGBMod, "bulk current model"),
        makeCelsiusParam("TNOM", &MosfetBSIM4Model::TNOM,
                         "reference temperature"),
        makeCelsiusParam("TEMP", &MosfetBSIM4Model::TEMP, "card temperature"),
    });
}

void MosfetBSIM4Model::applyGuards() {
    EG   = std::max(EG, 0.1);
    PHI  = std::max(PHI, 0.05);
    PB   = std::max(PB, 0.05);
    TOX  = std::max(TOX, 1e-10);
    U0   = std::max(U0, 1e-6);
    VSAT = std::max(VSAT, 1e3);
    UA   = std::max(UA, 0.0);
    UB   = std::max(UB, 0.0);
    GMIN = std::max(GMIN, 1e-15);
}

void MosfetBSIM4Model::onParamGiven(const std::string& key) {
    given_.insert(key);
}


void MosfetBSIM4Model::checkAndFixVersion() {
    // e-XSpice implements BSIM4 4.8.2. Version selection not supported.
    // If a VERSION param is provided, warn if not 4.82.
    // (No-op in our architecture; version is fixed at 4.8.2.)
}

bool MosfetBSIM4Model::processInstanceParams() {
    // Instance params are validated in Mosfet constructor.
    // This is a no-op in our architecture.
    return true;
}

void MosfetBSIM4Model::setupVersionPointers() {
    // Version-specific function pointers not needed (single version 4.8.2).
    // No-op.
}

void MosfetBSIM4Model::clearTemperatureData() {
    // Clear cached thermal data. In our architecture, thermal_ is recomputed
    // on demand via scaleTo(). No-op.
}

void MosfetBSIM4Model::printOutInstances() {
    // Debug: print instance info. No-op in production.
}
void MosfetBSIM4Model::applyDefaults() {
    // Xyce Model::processParams4p82_ (N_DEV_MOSFET_B4p82.C:7808-8106): junction
    // clamps first, then defaults computed from other parameters, applied only
    // when the card did not give them. Binning (L/W/P) and noise are not ported.
    const bool n = nmos_;
    if (PBS < 0.1) PBS = 0.1;
    if (PBSWS < 0.1) PBSWS = 0.1;
    if (PBSWGS < 0.1) PBSWGS = 0.1;
    if (PBD < 0.1) PBD = 0.1;
    if (PBSWD < 0.1) PBSWD = 0.1;
    if (PBSWGD < 0.1) PBSWGD = 0.1;

    if (!isGiven("TOX")) TOX = TOXE;
    if (!isGiven("TOXP")) TOXP = TOXE;
    if (!isGiven("TOXM")) TOXM = TOXE;
    if (!isGiven("DSUB")) DSUB = DROUT;
    if (!isGiven("VTH0")) VTH0 = n ? 0.7 : -0.7;
    if (!isGiven("VDDEOT")) VDDEOT = n ? 1.5 : -1.5;
    if (!isGiven("EU")) EU = n ? 1.67 : 1.0;
    if (!isGiven("UCS")) UCS = n ? 1.67 : 1.0;
    // VERSION 4.8.2 branch (> 4.80); MOBMOD selects the unit set.
    if (!isGiven("UA")) UA = (MOBMOD == 2 || MOBMOD == 6) ? 1.0e-15 : 1.0e-9;
    if (!isGiven("UC")) UC = (MOBMOD == 1 || MOBMOD == 5) ? -0.0465 : -0.0465e-9;
    if (!isGiven("UC1")) UC1 = (MOBMOD == 1 || MOBMOD == 5) ? -0.056 : -0.056e-9;
    if (!isGiven("U0")) { if (isGiven("KP")) U0 = KP / (EPSROX * 8.854214871e-12 / TOX); else U0 = n ? 0.067 : 0.025; }

    // 4.8.2 overrides the addPar default of FGIDL (0.0 -> 1.0).
    if (!isGiven("FGIDL")) FGIDL = 1.0;
    if (!isGiven("AGISL")) AGISL = AGIDL;
    if (!isGiven("BGISL")) BGISL = BGIDL;
    if (!isGiven("CGISL")) CGISL = CGIDL;
    if (!isGiven("EGISL")) EGISL = EGIDL;
    if (!isGiven("RGISL")) RGISL = RGIDL;
    if (!isGiven("KGISL")) KGISL = KGIDL;
    if (!isGiven("FGISL")) FGISL = FGIDL;

    if (!isGiven("AIGC")) AIGC = n ? 1.36e-2 : 9.80e-3;
    if (!isGiven("BIGC")) BIGC = n ? 1.71e-3 : 7.59e-4;
    if (!isGiven("CIGC")) CIGC = n ? 0.075 : 0.03;
    if (isGiven("AIGSD")) {
        AIGS = AIGD = AIGSD;
    } else {
        AIGSD = n ? 1.36e-2 : 9.80e-3;
        if (!isGiven("AIGS")) AIGS = AIGSD;
        if (!isGiven("AIGD")) AIGD = AIGSD;
    }
    if (isGiven("BIGSD")) {
        BIGS = BIGD = BIGSD;
    } else {
        BIGSD = n ? 1.71e-3 : 7.59e-4;
        if (!isGiven("BIGS")) BIGS = BIGSD;
        if (!isGiven("BIGD")) BIGD = BIGSD;
    }
    if (isGiven("CIGSD")) {
        CIGS = CIGD = CIGSD;
    } else {
        CIGSD = n ? 0.075 : 0.03;
        if (!isGiven("CIGS")) CIGS = CIGSD;
        if (!isGiven("CIGD")) CIGD = CIGSD;
    }

    if (!isGiven("IJTHDFWD")) IJTHDFWD = IJTHSFWD;
    if (!isGiven("IJTHDREV")) IJTHDREV = IJTHSREV;
    if (!isGiven("XJBVD")) XJBVD = XJBVS;
    if (!isGiven("BVD")) BVD = BVS;
    if (!isGiven("DMCI")) DMCI = DMCG;

    // Drain-side junction parameters follow the source side.
    if (!isGiven("JSD")) JSD = JSS;
    if (!isGiven("JSWD")) JSWD = JSWS;
    if (!isGiven("JSWGD")) JSWGD = JSWGS;
    if (!isGiven("PBD")) PBD = PBS;
    if (!isGiven("NJD")) NJD = NJS;
    if (!isGiven("XTID")) XTID = XTIS;
    if (!isGiven("MJD")) MJD = MJS;
    if (!isGiven("MJSWD")) MJSWD = MJSWS;
    if (!isGiven("MJSWGS")) MJSWGS = MJSWS;
    if (!isGiven("MJSWGD")) MJSWGD = MJSWGS;
    if (!isGiven("PBSWD")) PBSWD = PBSWS;
    if (!isGiven("PBSWGS")) PBSWGS = PBSWS;
    if (!isGiven("PBSWGD")) PBSWGD = PBSWD;
    if (!isGiven("CJD")) CJD = CJS;
    if (!isGiven("CJSWD")) CJSWD = CJSWS;
    if (!isGiven("CJSWGS")) CJSWGS = CJSWS;
    if (!isGiven("CJSWGD")) CJSWGD = CJSWGS;

    if (!isGiven("JTSD")) JTSD = JTSS;
    if (!isGiven("JTSSWD")) JTSSWD = JTSSWS;
    if (!isGiven("JTSSWGD")) JTSSWGD = JTSSWGS;
    if (!isGiven("NJTSD") && isGiven("NJTS")) NJTSD = NJTS;
    if (!isGiven("NJTSSWD") && isGiven("NJTSSW")) NJTSSWD = NJTSSW;
    if (!isGiven("NJTSSWGD") && isGiven("NJTSSWG")) NJTSSWGD = NJTSSWG;
    if (!isGiven("XTSD")) XTSD = XTSS;
    if (!isGiven("XTSSWD")) XTSSWD = XTSSWS;
    if (!isGiven("XTSSWGD")) XTSSWGD = XTSSWGS;
    if (!isGiven("TNJTSD") && isGiven("TNJTS")) TNJTSD = TNJTS;
    if (!isGiven("TNJTSSWD") && isGiven("TNJTSSW")) TNJTSSWD = TNJTSSW;
    if (!isGiven("TNJTSSWGD") && isGiven("TNJTSSWG")) TNJTSSWGD = TNJTSSWG;
    if (!isGiven("VTSD")) VTSD = VTSS;
    if (!isGiven("VTSSWD")) VTSSWD = VTSSWS;
    if (!isGiven("VTSSWGD")) VTSSWGD = VTSSWGS;
}

bool MosfetBSIM4Model::processParams(std::string& message) {
    if (LEVEL != 2) {
        message = "Only LEVEL=2 (BSIM4) is supported on this card";
        return false;
    }
    if (KP <= 0.0) { message = "KP must be positive"; return false; }
    if (TOX <= 0.0) { message = "TOX must be positive"; return false; }
    if (PHI <= 0.0) { message = "PHI must be positive"; return false; }
    if (TEMP <= 0.0 || TNOM <= 0.0) {
        message = "TEMP/TNOM must be positive (Kelvin)";
        return false;
    }
    checkAndFixVersion();
    applyDefaults();
    applyGuards();
    processInstanceParams();
    updateTemperature(TEMP);
    return true;
}

void MosfetBSIM4Model::updateTemperature(double tempKelvin) {
    thermal_ = scaleTo(tempKelvin);
}

MosfetBSIM4Thermal MosfetBSIM4Model::computeThermal(double tempKelvin) const {
    return scaleTo(tempKelvin);
}

MosfetBSIM4Thermal MosfetBSIM4Model::scaleTo(double Temp) const {
    MosfetBSIM4Thermal t;

    t.vt  = Const::KoverQ * Temp;
    double TRatio = Temp / TNOM;

    // Bandgap voltage (SPICE3 bandgap shift).
    double egfet = EG - (7.02e-4 * Temp * Temp) / (Temp + 1108.0);
    t.eg = egfet;

    // Intrinsic carrier concentration (ni).
    t.ni = 1.45e10 * (Temp / 300.15) * std::sqrt(Temp / 300.15)
         * std::exp(EG / (2.0 * Const::KoverQ) *
                    (1.0 / 300.15 - 1.0 / Temp));

    // Zero-bias mobility (temperature dependence via UTE).
    t.u0 = U0 * std::pow(TRatio, -UTE);

    // Saturation velocity (temperature dependence).
    t.vsat = VSAT * (1.0 - 2.0e-3 * (Temp - 300.15));

    // Body-effect coefficient scaled by TOX.
    t.k1ox = GAMMA * TOX / (Const::KoverQ * Temp);

    // Temperature coefficients for threshold (kept for compatibility).
    t.kt1 = 0.0;
    t.kt1l = 0.0;
    t.kt2 = 0.0;
    t.theta0vb0 = 0.0;

    // Mobility degradation coefficients.
    t.ua = UA + UA1 * (Temp - TNOM);
    t.ub = UB + UB1 * (Temp - TNOM);
    t.uc = UC + UC1 * (Temp - TNOM);
    t.ud = UD + UD1 * (Temp - TNOM);

    // Short-channel threshold coefficients.
    t.dvt0 = DVT0;
    t.dvt1 = DVT1;
    t.dvt2 = DVT2;
    t.dvt2w = DVT2W;

    t.gmin = GMIN;
    return t;
}

RCP<DeviceModel> MosfetBSIM4Model::fromCard(const std::string& name,
                                           const ParamMap& params,
                                           bool nmos) {
    ModelCard card;
    card.name = name;
    card.type = nmos ? "NMOS" : "PMOS";
    card.params = params;

    RCP<MosfetBSIM4Model> m = make_rcp<MosfetBSIM4Model>(name, nmos);
    m->LEVEL = static_cast<int>(card.get("LEVEL", 2.0));
    m->VTH0 = card.get("VTH0", m->VTH0);
    m->KP = card.get("KP", m->KP);
    m->LAMBDA = card.get("LAMBDA", m->LAMBDA);
    m->RDSW = card.get("RDSW", m->RDSW);
    m->PRWG = card.get("PRWG", m->PRWG);
    m->PRWB = card.get("PRWB", m->PRWB);
    m->PRT = card.get("PRT", m->PRT);
    m->RBDB = card.get("RBDB", m->RBDB);
    m->RBSB = card.get("RBSB", m->RBSB);
    m->RBPB = card.get("RBPB", m->RBPB);
    m->RBPS = card.get("RBPS", m->RBPS);
    m->RBPD = card.get("RBPD", m->RBPD);
    m->XRCRG1 = card.get("XRCRG1", m->XRCRG1);
    m->XRCRG2 = card.get("XRCRG2", m->XRCRG2);
    m->PHI = card.get("PHI", m->PHI);
    m->GAMMA = card.get("GAMMA", m->GAMMA);
    m->K2 = card.get("K2", m->K2);
    m->CBD = card.get("CBD", m->CBD);
    m->CBS = card.get("CBS", m->CBS);
    m->IS = card.get("IS", m->IS);
    m->PB = card.get("PB", m->PB);
    m->MJ = card.get("MJ", m->MJ);
    m->FC = card.get("FC", m->FC);
    m->CGSO = card.get("CGSO", m->CGSO);
    m->CGDO = card.get("CGDO", m->CGDO);
    m->CGBO = card.get("CGBO", m->CGBO);
    m->TOX = card.get("TOX", m->TOX);
    m->LD = card.get("LD", m->LD);
    m->EG = card.get("EG", m->EG);
    m->XTI = card.get("XTI", m->XTI);
    m->TNOM = card.get("TNOM", m->TNOM);
    m->TEMP = card.get("TEMP", m->TEMP);
    m->GMIN = card.get("GMIN", m->GMIN);
    // BSIM4-specific
    m->U0 = card.get("U0", m->U0);
    m->VSAT = card.get("VSAT", m->VSAT);
    m->UA = card.get("UA", m->UA);
    m->UB = card.get("UB", m->UB);
    m->UC = card.get("UC", m->UC);
    m->UD = card.get("UD", m->UD);
    m->A0 = card.get("A0", m->A0);
    m->B0 = card.get("B0", m->B0);
    m->B1 = card.get("B1", m->B1);
    m->KETA = card.get("KETA", m->KETA);
    m->KT1 = card.get("KT1", m->KT1);
    m->KT2 = card.get("KT2", m->KT2);
    m->ALPHA0 = card.get("ALPHA0", m->ALPHA0);
    m->AGIDL = card.get("AGIDL", m->AGIDL);
    m->BGIDL = card.get("BGIDL", m->BGIDL);
    m->CGIDL = card.get("CGIDL", m->CGIDL);
    m->EGIDL = card.get("EGIDL", m->EGIDL);
    m->AIGC = card.get("AIGC", m->AIGC);
    m->BIGC = card.get("BIGC", m->BIGC);
    m->CIGC = card.get("CIGC", m->CIGC);
    m->AIGS = card.get("AIGS", m->AIGS);
    m->AIGD = card.get("AIGD", m->AIGD);
    m->AIGBACC = card.get("AIGBACC", m->AIGBACC);
    m->AIGBINV = card.get("AIGBINV", m->AIGBINV);
    m->NIGC = card.get("NIGC", m->NIGC);
    m->TRNQSMOD = static_cast<int>(card.get("TRNQSMOD", m->TRNQSMOD));
    m->ACNQSMOD = static_cast<int>(card.get("ACNQSMOD", m->ACNQSMOD));
    m->TAUN = card.get("TAUN", m->TAUN);
    m->LMIN = card.get("LMIN", m->LMIN);
    m->LMAX = card.get("LMAX", m->LMAX);
    m->WMIN = card.get("WMIN", m->WMIN);
    m->WMAX = card.get("WMAX", m->WMAX);
    m->DMCG = card.get("DMCG", m->DMCG);
    m->DMCI = card.get("DMCI", m->DMCI);
    m->DMDG = card.get("DMDG", m->DMDG);
    m->SA = card.get("SA", m->SA);
    m->SB = card.get("SB", m->SB);
    m->SD = card.get("SD", m->SD);
    m->SCA = card.get("SCA", m->SCA);
    m->SCB = card.get("SCB", m->SCB);
    m->SCC = card.get("SCC", m->SCC);
    m->SC = card.get("SC", m->SC);
    m->KVTH0WE = card.get("KVTH0WE", m->KVTH0WE);
    m->K2WE = card.get("K2WE", m->K2WE);
    m->KU0WE = card.get("KU0WE", m->KU0WE);
    m->KF = card.get("KF", m->KF);
    m->NFACTOR = card.get("NFACTOR", m->NFACTOR);
    m->VOFF = card.get("VOFF", m->VOFF);
    m->AF = card.get("AF", m->AF);
    m->EF = card.get("EF", m->EF);
    m->RNOIA = card.get("RNOIA", m->RNOIA);
    m->RNOIB = card.get("RNOIB", m->RNOIB);
    m->RNOIC = card.get("RNOIC", m->RNOIC);
    m->AGISL = card.get("AGISL", m->AGISL);
    m->BGISL = card.get("BGISL", m->BGISL);
    m->CGISL = card.get("CGISL", m->CGISL);
    m->EGISL = card.get("EGISL", m->EGISL);
    m->ALPHA1 = card.get("ALPHA1", m->ALPHA1);
    m->UTE = card.get("UTE", m->UTE);
    m->DVT0 = card.get("DVT0", m->DVT0);
    m->DVT1 = card.get("DVT1", m->DVT1);
    m->DVT2 = card.get("DVT2", m->DVT2);
    m->DVT2W = card.get("DVT2W", m->DVT2W);
    m->DVTP0 = card.get("DVTP0", m->DVTP0);
    m->W0 = card.get("W0", m->W0);
    m->NFS = card.get("NFS", m->NFS);
    // Xyce BSIM4 4.8.2 addPar (N_DEV_MOSFET_B4.C), non-binned
    m->EOT = card.get("EOT", m->EOT);
    m->VDDEOT = card.get("VDDEOT", m->VDDEOT);
    m->TEMPEOT = card.get("TEMPEOT", m->TEMPEOT);
    m->LEFFEOT = card.get("LEFFEOT", m->LEFFEOT);
    m->WEFFEOT = card.get("WEFFEOT", m->WEFFEOT);
    m->ADOS = card.get("ADOS", m->ADOS);
    m->BDOS = card.get("BDOS", m->BDOS);
    m->TOXE = card.get("TOXE", m->TOXE);
    m->TOXP = card.get("TOXP", m->TOXP);
    m->TOXM = card.get("TOXM", m->TOXM);
    m->TOXREF = card.get("TOXREF", m->TOXREF);
    m->DTOX = card.get("DTOX", m->DTOX);
    m->EPSROX = card.get("EPSROX", m->EPSROX);
    m->CDSC = card.get("CDSC", m->CDSC);
    m->CDSCB = card.get("CDSCB", m->CDSCB);
    m->CDSCD = card.get("CDSCD", m->CDSCD);
    m->CIT = card.get("CIT", m->CIT);
    m->XJ = card.get("XJ", m->XJ);
    m->AT = card.get("AT", m->AT);
    m->AGS = card.get("AGS", m->AGS);
    m->A1 = card.get("A1", m->A1);
    m->A2 = card.get("A2", m->A2);
    m->PHIG = card.get("PHIG", m->PHIG);
    m->EPSRGATE = card.get("EPSRGATE", m->EPSRGATE);
    m->EASUB = card.get("EASUB", m->EASUB);
    m->EPSRSUB = card.get("EPSRSUB", m->EPSRSUB);
    m->NI0SUB = card.get("NI0SUB", m->NI0SUB);
    m->BG0SUB = card.get("BG0SUB", m->BG0SUB);
    m->TBGASUB = card.get("TBGASUB", m->TBGASUB);
    m->TBGBSUB = card.get("TBGBSUB", m->TBGBSUB);
    m->NSUB = card.get("NSUB", m->NSUB);
    m->NDEP = card.get("NDEP", m->NDEP);
    m->NSD = card.get("NSD", m->NSD);
    m->PHIN = card.get("PHIN", m->PHIN);
    m->NGATE = card.get("NGATE", m->NGATE);
    m->GAMMA1 = card.get("GAMMA1", m->GAMMA1);
    m->GAMMA2 = card.get("GAMMA2", m->GAMMA2);
    m->VBX = card.get("VBX", m->VBX);
    m->VBM = card.get("VBM", m->VBM);
    m->XT = card.get("XT", m->XT);
    m->K1 = card.get("K1", m->K1);
    m->KT1L = card.get("KT1L", m->KT1L);
    m->K3 = card.get("K3", m->K3);
    m->K3B = card.get("K3B", m->K3B);
    m->DVTP1 = card.get("DVTP1", m->DVTP1);
    m->DVTP2 = card.get("DVTP2", m->DVTP2);
    m->DVTP3 = card.get("DVTP3", m->DVTP3);
    m->DVTP4 = card.get("DVTP4", m->DVTP4);
    m->DVTP5 = card.get("DVTP5", m->DVTP5);
    m->LPE0 = card.get("LPE0", m->LPE0);
    m->LPEB = card.get("LPEB", m->LPEB);
    m->DVT0W = card.get("DVT0W", m->DVT0W);
    m->DVT1W = card.get("DVT1W", m->DVT1W);
    m->DROUT = card.get("DROUT", m->DROUT);
    m->DSUB = card.get("DSUB", m->DSUB);
    m->UA1 = card.get("UA1", m->UA1);
    m->UB1 = card.get("UB1", m->UB1);
    m->UC1 = card.get("UC1", m->UC1);
    m->UD1 = card.get("UD1", m->UD1);
    m->UP = card.get("UP", m->UP);
    m->LP = card.get("LP", m->LP);
    m->EU = card.get("EU", m->EU);
    m->UCS = card.get("UCS", m->UCS);
    m->UCSTE = card.get("UCSTE", m->UCSTE);
    m->MINV = card.get("MINV", m->MINV);
    m->VOFFL = card.get("VOFFL", m->VOFFL);
    m->XPART = card.get("XPART", m->XPART);
    m->DELTA = card.get("DELTA", m->DELTA);
    m->RSH = card.get("RSH", m->RSH);
    m->RDSWMIN = card.get("RDSWMIN", m->RDSWMIN);
    m->RSW = card.get("RSW", m->RSW);
    m->RDW = card.get("RDW", m->RDW);
    m->RDWMIN = card.get("RDWMIN", m->RDWMIN);
    m->RSWMIN = card.get("RSWMIN", m->RSWMIN);
    m->ETA0 = card.get("ETA0", m->ETA0);
    m->ETAB = card.get("ETAB", m->ETAB);
    m->PCLM = card.get("PCLM", m->PCLM);
    m->PDIBLC1 = card.get("PDIBLC1", m->PDIBLC1);
    m->PDIBLC2 = card.get("PDIBLC2", m->PDIBLC2);
    m->PDIBLCB = card.get("PDIBLCB", m->PDIBLCB);
    m->FPROUT = card.get("FPROUT", m->FPROUT);
    m->PDITS = card.get("PDITS", m->PDITS);
    m->PDITSL = card.get("PDITSL", m->PDITSL);
    m->PDITSD = card.get("PDITSD", m->PDITSD);
    m->PSCBE1 = card.get("PSCBE1", m->PSCBE1);
    m->PSCBE2 = card.get("PSCBE2", m->PSCBE2);
    m->PVAG = card.get("PVAG", m->PVAG);
    m->JSS = card.get("JSS", m->JSS);
    m->JSWS = card.get("JSWS", m->JSWS);
    m->JSWGS = card.get("JSWGS", m->JSWGS);
    m->PBS = card.get("PBS", m->PBS);
    m->NJS = card.get("NJS", m->NJS);
    m->XTIS = card.get("XTIS", m->XTIS);
    m->MJS = card.get("MJS", m->MJS);
    m->PBSWS = card.get("PBSWS", m->PBSWS);
    m->MJSWS = card.get("MJSWS", m->MJSWS);
    m->PBSWGS = card.get("PBSWGS", m->PBSWGS);
    m->MJSWGS = card.get("MJSWGS", m->MJSWGS);
    m->CJS = card.get("CJS", m->CJS);
    m->CJSWS = card.get("CJSWS", m->CJSWS);
    m->CJSWGS = card.get("CJSWGS", m->CJSWGS);
    m->JSD = card.get("JSD", m->JSD);
    m->JSWD = card.get("JSWD", m->JSWD);
    m->JSWGD = card.get("JSWGD", m->JSWGD);
    m->PBD = card.get("PBD", m->PBD);
    m->NJD = card.get("NJD", m->NJD);
    m->XTID = card.get("XTID", m->XTID);
    m->MJD = card.get("MJD", m->MJD);
    m->PBSWD = card.get("PBSWD", m->PBSWD);
    m->MJSWD = card.get("MJSWD", m->MJSWD);
    m->PBSWGD = card.get("PBSWGD", m->PBSWGD);
    m->MJSWGD = card.get("MJSWGD", m->MJSWGD);
    m->CJD = card.get("CJD", m->CJD);
    m->CJSWD = card.get("CJSWD", m->CJSWD);
    m->CJSWGD = card.get("CJSWGD", m->CJSWGD);
    m->TPB = card.get("TPB", m->TPB);
    m->TCJ = card.get("TCJ", m->TCJ);
    m->TPBSW = card.get("TPBSW", m->TPBSW);
    m->TCJSW = card.get("TCJSW", m->TCJSW);
    m->TPBSWG = card.get("TPBSWG", m->TPBSWG);
    m->TCJSWG = card.get("TCJSWG", m->TCJSWG);
    m->MOIN = card.get("MOIN", m->MOIN);
    m->NOFF = card.get("NOFF", m->NOFF);
    m->XGW = card.get("XGW", m->XGW);
    m->XGL = card.get("XGL", m->XGL);
    m->RSHG = card.get("RSHG", m->RSHG);
    m->NGCON = card.get("NGCON", m->NGCON);
    m->VTL = card.get("VTL", m->VTL);
    m->LC = card.get("LC", m->LC);
    m->XN = card.get("XN", m->XN);
    m->VFBSDOFF = card.get("VFBSDOFF", m->VFBSDOFF);
    m->TVFBSDOFF = card.get("TVFBSDOFF", m->TVFBSDOFF);
    m->TVOFF = card.get("TVOFF", m->TVOFF);
    m->TNFACTOR = card.get("TNFACTOR", m->TNFACTOR);
    m->TETA0 = card.get("TETA0", m->TETA0);
    m->LINT = card.get("LINT", m->LINT);
    m->LL = card.get("LL", m->LL);
    m->LLN = card.get("LLN", m->LLN);
    m->LW = card.get("LW", m->LW);
    m->LWN = card.get("LWN", m->LWN);
    m->LWL = card.get("LWL", m->LWL);
    m->WR = card.get("WR", m->WR);
    m->WINT = card.get("WINT", m->WINT);
    m->DWG = card.get("DWG", m->DWG);
    m->DWB = card.get("DWB", m->DWB);
    m->WL = card.get("WL", m->WL);
    m->WLN = card.get("WLN", m->WLN);
    m->WW = card.get("WW", m->WW);
    m->WWN = card.get("WWN", m->WWN);
    m->WWL = card.get("WWL", m->WWL);
    m->XW = card.get("XW", m->XW);
    m->XL = card.get("XL", m->XL);
    m->BETA0 = card.get("BETA0", m->BETA0);
    m->RGIDL = card.get("RGIDL", m->RGIDL);
    m->KGIDL = card.get("KGIDL", m->KGIDL);
    m->FGIDL = card.get("FGIDL", m->FGIDL);
    m->RGISL = card.get("RGISL", m->RGISL);
    m->KGISL = card.get("KGISL", m->KGISL);
    m->FGISL = card.get("FGISL", m->FGISL);
    m->AIGSD = card.get("AIGSD", m->AIGSD);
    m->BIGSD = card.get("BIGSD", m->BIGSD);
    m->CIGSD = card.get("CIGSD", m->CIGSD);
    m->BIGS = card.get("BIGS", m->BIGS);
    m->CIGS = card.get("CIGS", m->CIGS);
    m->BIGD = card.get("BIGD", m->BIGD);
    m->CIGD = card.get("CIGD", m->CIGD);
    m->BIGBACC = card.get("BIGBACC", m->BIGBACC);
    m->CIGBACC = card.get("CIGBACC", m->CIGBACC);
    m->BIGBINV = card.get("BIGBINV", m->BIGBINV);
    m->CIGBINV = card.get("CIGBINV", m->CIGBINV);
    m->NIGBINV = card.get("NIGBINV", m->NIGBINV);
    m->NIGBACC = card.get("NIGBACC", m->NIGBACC);
    m->NTOX = card.get("NTOX", m->NTOX);
    m->EIGBINV = card.get("EIGBINV", m->EIGBINV);
    m->PIGCD = card.get("PIGCD", m->PIGCD);
    m->POXEDGE = card.get("POXEDGE", m->POXEDGE);
    m->IJTHDFWD = card.get("IJTHDFWD", m->IJTHDFWD);
    m->IJTHSFWD = card.get("IJTHSFWD", m->IJTHSFWD);
    m->IJTHDREV = card.get("IJTHDREV", m->IJTHDREV);
    m->IJTHSREV = card.get("IJTHSREV", m->IJTHSREV);
    m->XJBVD = card.get("XJBVD", m->XJBVD);
    m->XJBVS = card.get("XJBVS", m->XJBVS);
    m->BVD = card.get("BVD", m->BVD);
    m->BVS = card.get("BVS", m->BVS);
    m->JTSS = card.get("JTSS", m->JTSS);
    m->JTSD = card.get("JTSD", m->JTSD);
    m->JTSSWS = card.get("JTSSWS", m->JTSSWS);
    m->JTSSWD = card.get("JTSSWD", m->JTSSWD);
    m->JTSSWGS = card.get("JTSSWGS", m->JTSSWGS);
    m->JTSSWGD = card.get("JTSSWGD", m->JTSSWGD);
    m->JTWEFF = card.get("JTWEFF", m->JTWEFF);
    m->NJTS = card.get("NJTS", m->NJTS);
    m->NJTSSW = card.get("NJTSSW", m->NJTSSW);
    m->NJTSSWG = card.get("NJTSSWG", m->NJTSSWG);
    m->NJTSD = card.get("NJTSD", m->NJTSD);
    m->NJTSSWD = card.get("NJTSSWD", m->NJTSSWD);
    m->NJTSSWGD = card.get("NJTSSWGD", m->NJTSSWGD);
    m->XTSS = card.get("XTSS", m->XTSS);
    m->XTSD = card.get("XTSD", m->XTSD);
    m->XTSSWS = card.get("XTSSWS", m->XTSSWS);
    m->XTSSWD = card.get("XTSSWD", m->XTSSWD);
    m->XTSSWGS = card.get("XTSSWGS", m->XTSSWGS);
    m->XTSSWGD = card.get("XTSSWGD", m->XTSSWGD);
    m->TNJTS = card.get("TNJTS", m->TNJTS);
    m->TNJTSSW = card.get("TNJTSSW", m->TNJTSSW);
    m->TNJTSSWG = card.get("TNJTSSWG", m->TNJTSSWG);
    m->TNJTSD = card.get("TNJTSD", m->TNJTSD);
    m->TNJTSSWD = card.get("TNJTSSWD", m->TNJTSSWD);
    m->TNJTSSWGD = card.get("TNJTSSWGD", m->TNJTSSWGD);
    m->VTSS = card.get("VTSS", m->VTSS);
    m->VTSD = card.get("VTSD", m->VTSD);
    m->VTSSWS = card.get("VTSSWS", m->VTSSWS);
    m->VTSSWD = card.get("VTSSWD", m->VTSSWD);
    m->VTSSWGS = card.get("VTSSWGS", m->VTSSWGS);
    m->VTSSWGD = card.get("VTSSWGD", m->VTSSWGD);
    m->GBMIN = card.get("GBMIN", m->GBMIN);
    m->RBPS0 = card.get("RBPS0", m->RBPS0);
    m->RBPSL = card.get("RBPSL", m->RBPSL);
    m->RBPSW = card.get("RBPSW", m->RBPSW);
    m->RBPSNF = card.get("RBPSNF", m->RBPSNF);
    m->RBPD0 = card.get("RBPD0", m->RBPD0);
    m->RBPDL = card.get("RBPDL", m->RBPDL);
    m->RBPDW = card.get("RBPDW", m->RBPDW);
    m->RBPDNF = card.get("RBPDNF", m->RBPDNF);
    m->RBPBX0 = card.get("RBPBX0", m->RBPBX0);
    m->RBPBXL = card.get("RBPBXL", m->RBPBXL);
    m->RBPBXW = card.get("RBPBXW", m->RBPBXW);
    m->RBPBXNF = card.get("RBPBXNF", m->RBPBXNF);
    m->RBPBY0 = card.get("RBPBY0", m->RBPBY0);
    m->RBPBYL = card.get("RBPBYL", m->RBPBYL);
    m->RBPBYW = card.get("RBPBYW", m->RBPBYW);
    m->RBPBYNF = card.get("RBPBYNF", m->RBPBYNF);
    m->RBSBX0 = card.get("RBSBX0", m->RBSBX0);
    m->RBSBY0 = card.get("RBSBY0", m->RBSBY0);
    m->RBDBX0 = card.get("RBDBX0", m->RBDBX0);
    m->RBDBY0 = card.get("RBDBY0", m->RBDBY0);
    m->RBSDBXL = card.get("RBSDBXL", m->RBSDBXL);
    m->RBSDBXW = card.get("RBSDBXW", m->RBSDBXW);
    m->RBSDBXNF = card.get("RBSDBXNF", m->RBSDBXNF);
    m->RBSDBYL = card.get("RBSDBYL", m->RBSDBYL);
    m->RBSDBYW = card.get("RBSDBYW", m->RBSDBYW);
    m->RBSDBYNF = card.get("RBSDBYNF", m->RBSDBYNF);
    m->SAREF = card.get("SAREF", m->SAREF);
    m->SBREF = card.get("SBREF", m->SBREF);
    m->WLOD = card.get("WLOD", m->WLOD);
    m->KU0 = card.get("KU0", m->KU0);
    m->KVSAT = card.get("KVSAT", m->KVSAT);
    m->KVTH0 = card.get("KVTH0", m->KVTH0);
    m->TKU0 = card.get("TKU0", m->TKU0);
    m->LLODKU0 = card.get("LLODKU0", m->LLODKU0);
    m->WLODKU0 = card.get("WLODKU0", m->WLODKU0);
    m->LLODVTH = card.get("LLODVTH", m->LLODVTH);
    m->WLODVTH = card.get("WLODVTH", m->WLODVTH);
    m->STK2 = card.get("STK2", m->STK2);
    m->LODK2 = card.get("LODK2", m->LODK2);
    m->STETA0 = card.get("STETA0", m->STETA0);
    m->LODETA0 = card.get("LODETA0", m->LODETA0);
    m->WEB = card.get("WEB", m->WEB);
    m->WEC = card.get("WEC", m->WEC);
    m->SCREF = card.get("SCREF", m->SCREF);
    m->GIDLCLAMP = card.get("GIDLCLAMP", m->GIDLCLAMP);
    m->IDOVVDS = card.get("IDOVVDS", m->IDOVVDS);
    m->CVCHARGEMOD = static_cast<int>(card.get("CVCHARGEMOD", m->CVCHARGEMOD));
    m->CAPMOD = static_cast<int>(card.get("CAPMOD", m->CAPMOD));
    m->DIOMOD = static_cast<int>(card.get("DIOMOD", m->DIOMOD));
    m->RDSMOD = static_cast<int>(card.get("RDSMOD", m->RDSMOD));
    m->MOBMOD = static_cast<int>(card.get("MOBMOD", m->MOBMOD));
    m->PERMOD = static_cast<int>(card.get("PERMOD", m->PERMOD));
    m->GEOMOD = static_cast<int>(card.get("GEOMOD", m->GEOMOD));
    m->RGEOMOD = static_cast<int>(card.get("RGEOMOD", m->RGEOMOD));
    m->FNOIMOD = static_cast<int>(card.get("FNOIMOD", m->FNOIMOD));
    m->TNOIMOD = static_cast<int>(card.get("TNOIMOD", m->TNOIMOD));
    m->MTRLMOD = static_cast<int>(card.get("MTRLMOD", m->MTRLMOD));
    m->MTRLCOMPATMOD = static_cast<int>(card.get("MTRLCOMPATMOD", m->MTRLCOMPATMOD));
    m->TEMPMOD = static_cast<int>(card.get("TEMPMOD", m->TEMPMOD));
    m->GIDLMOD = static_cast<int>(card.get("GIDLMOD", m->GIDLMOD));
    m->PARAMCHK = static_cast<int>(card.get("PARAMCHK", m->PARAMCHK));
    m->BINUNIT = static_cast<int>(card.get("BINUNIT", m->BINUNIT));
    m->DMCGT = card.get("DMCGT", m->DMCGT);
    m->RBODYMOD = static_cast<int>(card.get("RBODYMOD", m->RBODYMOD));
    m->RGATEMOD = static_cast<int>(card.get("RGATEMOD", m->RGATEMOD));
    m->WPEMOD = static_cast<int>(card.get("WPEMOD", m->WPEMOD));
    m->IGCMod = static_cast<int>(card.get("IGCMOD", m->IGCMod));
    m->IGBMod = static_cast<int>(card.get("IGBMOD", m->IGBMod));
    for (const auto& kv : params) {
        if (m->hasParam(kv.first)) m->onParamGiven(kv.first);
    }
    if (card.given("TNOM")) m->TNOM = card.get("TNOM", 27.0) + Const::CtoK;
    if (card.given("TEMP")) m->TEMP = card.get("TEMP", 27.0) + Const::CtoK;
    return m;
}

RCP<const MosfetBSIM4Model> MosfetBSIM4Model::defaultNmos() {
    static RCP<MosfetBSIM4Model> instance;
    if (!instance) {
        instance = make_rcp<MosfetBSIM4Model>("$DEFAULT_NMOS_BSIM4", true);
        std::string message;
        if (!instance->processParams(message)) {
            throw std::logic_error("default NMOS BSIM4 model invalid: " +
                                   message);
        }
    }
    return instance;
}

RCP<const MosfetBSIM4Model> MosfetBSIM4Model::defaultPmos() {
    static RCP<MosfetBSIM4Model> instance;
    if (!instance) {
        instance = make_rcp<MosfetBSIM4Model>("$DEFAULT_PMOS_BSIM4", false);
        std::string message;
        if (!instance->processParams(message)) {
            throw std::logic_error("default PMOS BSIM4 model invalid: " +
                                   message);
        }
    }
    return instance;
}

} // namespace e_XSpice
