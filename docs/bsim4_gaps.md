# Verificação de lacunas: BSIM4 do Xyce vs. e-XSpice

Comparação entre `external_codes/Xyce/OpenModels/N_DEV_MOSFET_B4{,p82}.C` (BSIM4 4.8.2/4.8.3) e `src/models/core_models/MosfetBSIM4{Device,Model}.{cpp,h}` + `src/devices/Mosfet.{cpp,h}`. O que já foi feito no P0 (M/NF/área, U0, junções, nós internos, VFB, K1/K2) não é listado como lacuna.

## 1. Geometria / PAeffGeo / NumFingerDiff / RdseffGeo
**Status: parcial → melhorado no P1.**
- `NumFingerDiff` portado para `min==0` e `min==1`.
- `PAeffGeo` portado para todos os casos 0..10 (switch completo).
- `PERMOD` parcialmente (strip do gate-edge).
- `RdseffGeo`/`RdsEndIso`/`RdsEndSha` ainda não portados (aproximação linear).
- `GEOMOD` instância não honrado (usa model_->GEOMOD).

## 2. Temperatura (updateTemperature4p82)
**Status: parcial.**
- `scaleTo()` cobre UTE, junction sat-current (JSS/JSD/JSWS/...), C-V (CJS/CJD/...), `vcrit`.
- Faltam: `SjctTempExponent/SjctEmissionCoeff` usados no Xyce; `KT1/KT2` com `KT1L`; `AT` (vsat temp); `PRT` para Rds; `UCSTE`; `TVFBSDOFF`; `TVOFF/TNFACTOR/TETA0`; `TBGASUB/TBGBSUB` (bandgap) — e-XSpice usa `egfet` simplificado.
- `TEMPMOD` (self-heating) não portado: `tempMod`, `SHMOD`, `RTH0`, `CTH0`, `WTH0` ausentes; nenhuma equação de temperatura dinâmica.

## 3. Mobilidade
**Status: parcial.**
- `ueff` portado (universal/Coulombic/high-k, LOD, WPE).
- Faltam: `UA1/UB1/UC1` (mobilidade de 2ª ordem com Vds); `UP/LP` (mobilidade de canal curto); `EU` exponente; `UCS/UCSTE`; termo `RDSW` na mobilidade (bias-dependent dentro do device já usa RDSW/PRWG/PRWB, mas não a versão completa com `RSW/RDW` mínimos).

## 4. Vth / body effect / DIBL
**Status: parcial → melhorado no P1.**
- Vth portado com K1/K2, K3/K3B, KETA, DIBL (eta0/etab), Voff, NFACTOR, MINV.
- DIBL melhorado com DSUB/DROUT/FPROUT.
- CLM/PVAG portados (PCLM/PVAG).
- Faltam: `DVT0/DVT1/DVT2/DVT0W/DVT1W/DVT2W` completo; `DVTP0..5`; `LPE0/LPEB`; `VTH0` L/W dep; `CDSC/CDSCB/CDSCD/CIT`; `A1/A2`; `AGS`; `PDITS/PDITSD/PDITSL`; `DELTA`.

## 5. Junctions (saturation + C-V)
**Status: parcial → melhorado no P1.**
- Saturation: bottom + sidewall + gate-edge, temp-scaled, M-aplicado.
- Junction tunneling (JTSS/JTSD) portado (simplificado).
- Breakdown (BVD/BVS) portado (simplificado).
- `IJTHDFWD/IJTHSFWD/IJTHDREV/IJTHSREV` + `DioIjthVjmEval` ainda não portados.
- `GIDLCLAMP` portado.
- `NJTS` exponentes por termo, `XTSS/...`, `VTSS/...` parcialmente.
- C-V: intrinsic charges (capMod 0/1/2) portados parcialmente; `CVCHARGEMOD` (0/1/2) honrado parcialmente; `CAPMOD` (0/1/2/3) não totalmente.

## 6. Internal nodes (RDSMOD/RGATEMOD/RBODYMOD/RGEOMOD/GEOMOD/MIN)
**Status: parcial → melhorado no P1.**
- `RDSMOD` (bias-dependent, dentro do device) + sheet R (RSH*NRD/NRS) portados.
- `RGATEMOD==1` portado; `==2/3` desabilitados com warning.
- `RBODYMOD==1` portado (rede de 5 resistores com 1 nó interno: RBDB/RBSB/RBPB/RBPS/RBPD). `==2` cai para 1 com warning.
- `RGEOMOD` não portado.
- `GEOMOD` (minimize source/drain) honrado via MIN (instância).

## 7. C-V / cargas intrínsecas
**Status: parcial → melhorado no P1.**
- `ctmCharges` (capMod 0) e `capMod` 1/2 portados.
- `POLYDEPLETION` (`polyDepletion`) portado.
- `CLC/CLE` (capacitância de franja) portado.
- `DLC/DWC` (overlap) portado.
- `CKAPPAS/CKAPPAD/CF` portados.
- `VFB`/`VFBCV`/`VOFFCV`/`VOFFCVL`/`MINVCV`/`ACDE` parcialmente.

## 8. NQS (TRNQSMOD/ACNQSMOD)
**Status: não portado.**
- Parâmetros declarados mas nenhuma equação. Faltam: `TAUN`, `ELM`, `XRCRG1/XRCRG2` (gate resistance distribuído), `TRNQSMOD` (0/1), `ACNQSMOD` (0/1), `QDEF` (carga de defeito), e as equações de relaxação.

## 9. Noise
**Status: parcial.**
- `eval1ovFNoise` declarado mas **não definido** (link latente). `getNumNoiseSources` declarado mas não definido.
- Faltam: `FNOIMOD/TNOIMOD`, `KF/AF/EF`, `RNOIA/RNOIB/RNOIC`, `TNOIA/TNOIB/TNOIC`, `NTNOI`, `LINTNOI`, e o modelo de ruído térmico completo (canal + gate + body).

## 10. Gate tunneling (IGC/IGB) e GIDL/GISL
**Status: parcial → melhorado no P1.**
- `AIGC/BIGC/CIGC/AIGS/...` portados; `CIGC` portado.
- `GIDLCLAMP` portado.
- `GIDL` com `AGIDL/BGIDL/CGIDL/EGIDL/FGIDL/KGIDL/RGIDL` portado.
- `GISL` com `AGISL/BGISL/CGISL/EGISL/FGISL/KGISL/RGISL` portado.
- `IDOVVDSC` ainda não portado.

## 11. Body diode / breakdown
**Status: parcial.**
- Diode básico (saturation) portado.
- Faltam: `BV`/`IBV` (breakdown), `DIO` model variation, `XJBVD/XJBVS`, `IJTH*`, `VTSS/VTSD/...` (por termo), `NTNOI`.

## 12. DIBL/CLM/velocity saturation
**Status: parcial.**
- `PCLM/PDIBL1/PDIBL2` não portados; `DROUT/DSUB` não portados; `A1/A2` não portados; `PVAG` portado parcialmente.

## 13. High-k / metal gate / poly depletion
**Status: não portado.**
- `PHIG/EPSRGATE/EASUB/EPSRSUB/NI0SUB/BG0SUB/TBGASUB/TBGBSUB/NDPE/NSUB/NSD` parcialmente.
- `POLYDEPLETION` não portado.
- `EOT/VDDEOT/TEMPEOT/LEFFEOT/WEFFEOT/ADOS/BDOS` (EOT scaling) não portados.

## 14. Binagem / corners
**Status: não portado.**
- `LMIN/LMAX/WMIN/WMAX` e binagem de parâmetros (L/W/NF) não portados.

## 15. Outros parâmetros faltantes
- `TNJTS/TNJTSWS/...` (temp exponentes por termo de junção) — parcialmente via `XTIS/XTID`.
- `VTSS/VTSD/VTSSWS/...` (junction potentials por termo) — parcialmente.
- `JTWEFF`, `CKAPPAS/CKAPPAD/CF` (cap de franja/overlap).
- `TVFBSDOFF`, `TVFBOFF`, `TVOFFCV`, `MINVCV`.
- `GMIN` (paralelo) portado.

## Resumo quantitativo
- Parâmetros do model Xyce: ~264 membros principais. e-XSpice declara ~432 (inclui ruído/AC/overlap extras). Faltam ~21 parâmetros específicos (tnoia/tnoib/tnoic/ntnoi/lintnoi/gidlclamp/idovvdsc/vfbsdoff/pdibl1/2/b/tvoffcv/modType/versionDouble).
- Funções do Xyce Instance: ~30 (processParams, updateTemperature, updateIntermediateVars, updatePrimaryState, auxChargeCalculations, setupCapacitors, loadDAE*, PAeffGeo, NumFingerDiff, RdseffGeo, RdsEndIso, RdsEndSha, polyDepletion, DioIjthVjmEval, Eval1ovFNoise, getNoiseSources, register*). e-XSpice implementa um subconjunto (evaluateBase, loadF/Q/Jacobian, cvParams, auxCharge parcial).

## P2 (implementado nesta entrega)
- `eval1ovFNoise` e `getNumNoiseSources` definidos (link latente resolvido).
- Noise completo: thermal (RNOIA/RNOIB/RNOIC), gate/source thermal, flicker full (KF/AF/EF e NOIA/NOIB/NOIC).
- Self-heating (TEMPMOD): equilíbrio térmico estático T = Tamb + RTH0*P com re-avaliação iterativa.
- `TRNQSMOD`/`ACNQSMOD` desabilitados com warning (NQS não portado — requer infraestrutura de estado DAE não disponível).
- Binagem LMIN/LMAX/WMIN/WMAX honrada no construtor do Mosfet.
- High-k/metal gate params declarados; VFB usado no C-V.

## Prioridade sugerida (P1/P2)
1. **P1 (física core):** DIBL completo (DVT*, DSUB, DROUT), CLM (A1/A2, PCLM, PDIBL*), velocity saturation completo, poly depletion, C-V completo (capMod 0/1/2/3, overlap DLC/DWC, franja), junction tunneling + breakdown, body diode completo.
2. **P1 (geometria):** PAeffGeo casos 1/2/3, RdseffGeo completo, RdsEndIso/Sha, MIN, GEOMOD.
3. **P2 (avançado):** NQS (TRNQS/ACNQS), self-heating (TEMPMOD), noise completo, gate tunneling completo, high-k/metal gate, binagem.
4. **P2 (infra):** definir `eval1ovFNoise` e `getNumNoiseSources` (link latente), honrar `TRNQSMOD/ACNQSMOD` ou desabilitar com warning.
