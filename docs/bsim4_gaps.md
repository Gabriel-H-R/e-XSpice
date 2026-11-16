# Verificação de lacunas: BSIM4 do Xyce vs. e-XSpice

Comparação sistemática entre `external_codes/Xyce/OpenModels/N_DEV_MOSFET_B4{,p82}.C` (BSIM4 4.8.2/4.8.3) e `src/models/core_models/MosfetBSIM4{Device,Model}.{cpp,h}` + `src/devices/Mosfet.{cpp,h}`. Este documento é o resultado da verificação solicitada: ele lista, com referências de função/linha do Xyce, o que **falta ou está simplificado** no modelo físico do e-XSpice.

Legenda de status:
- **OK** – portado com fidelidade Xyce.
- **PARCIAL** – portado mas simplificado / com termos faltantes.
- **NÃO** – não portado (parâmetro declarado mas sem equação, ou ausente).
- **WARN** – modo desabilitado com aviso (fallback seguro).

---

## 1. Cobertura de funções do Xyce Instance

| Função Xyce (N_DEV_MOSFET_B4p82.C) | e-XSpice | Status |
|---|---|---|
| `Instance::processParams4p82_` (7808) | `MosfetBSIM4Model::applyDefaults()` + `processParams()` | PARCIAL |
| `Instance::updateTemperature4p82_` (346) | `MosfetBSIM4Model::scaleTo()` | PARCIAL |
| `Instance::updateIntermediateVars4p82_` (2744) | `MosfetBSIM4Device::evaluateBase()` | PARCIAL |
| `Instance::updatePrimaryState` (6881) | charges inline em `evaluateBase` | OK (equivalente) |
| `Instance::auxChargeCalculations` (7290) | C-V inline | OK |
| `Instance::setupCapacitors_newDAE` (7412) | C-V inline | OK |
| `Instance::loadDAEFVector` (7899) | `Mosfet::loadF` + `MosfetBSIM4Device::loadF` | OK |
| `Instance::loadDAEQVector` (7125) | `Mosfet::loadQ` + device `loadQ` | OK |
| `Instance::loadDAEdFdx` (8541) | `loadJacobianF` | OK |
| `Instance::loadDAEdQdx` (8431) | `loadJacobianQ` (diferenças centrais) | OK |
| `Instance::polyDepletion` (8819) | `MosfetBSIM4Device::polyDepletion` | OK (mas não chamado no fluxo principal – NGATE não portado) |
| `Instance::NumFingerDiff` (8889) | `MosfetBSIM4Device::numFingerDiff` | OK |
| `Instance::PAeffGeo` (8931) | `MosfetBSIM4Device::paEffGeo` | OK (todos os casos 0..10) |
| `Instance::RdseffGeo` | **não implementado** | NÃO |
| `Instance::RdsEndIso` / `RdsEndIso4p82_` | **não implementado** (aproximação linear) | NÃO |
| `Instance::RdsEndSha` | **não implementado** | NÃO |
| `Instance::DioIjthVjmEval` | **não implementado** | NÃO |
| `Instance::Eval1ovFNoise` (B4.C) | `MosfetBSIM4Device::eval1ovFNoise` | PARCIAL |
| `Instance::getNumNoiseSources` | `MosfetBSIM4Device::getNumNoiseSources` | PARCIAL (retorna 2 fixo) |
| `Instance::getNoiseSources4p82_` | `Mosfet::loadNoiseSources` | PARCIAL |
| `Instance::setupNoiseSources4p82_` | não existe (estático) | NÃO |

**Funções do Master/Model do Xyce (não portadas por serem infraestrutura de solver):** `Master::updateState`, `Master::loadDAEVectors`, `Master::loadDAEMatrices`, `Model::clearTemperatureData`, `Model::forEachInstance`. São específicas do solver Xyce e não aplicáveis à arquitetura do e-XSpice.

---

## 2. Blocos físicos de `updateIntermediateVars4p82_` (Xyce) vs. `evaluateBase` (e-XSpice)

| Bloco Xyce (linha aprox.) | e-XSpice | Status / o que falta |
|---|---|---|
| Junction diode DC model (3348-3695) | sat. corrente + breakdown + tunneling | PARCIAL – faltam `IJTHDFWD/SFWD/DREV/SREV`, `XJBVD/XJBVS`, `DioIjthVjmEval`, `VTSS/VTSD/...` por termo, `NJTS` exponents por termo (parcial via XTIS/XTID) |
| Vth – narrow width (3795) | Vth_NarrowW = TOXE*PHI/(weff+W0) | OK |
| Vth – body effect K1/K2 (3894) | t.k1ox/t.k2ox | OK |
| Vth – K3/K3B (3795) | model_->K3/K3B | OK |
| Vth – KETA (3810) | model_->KETA | OK |
| Vth – DVT0/DVT1/DVT2 (3810) | simplificado (DVT0/DVT1/DVT2, DVT0W ausente) | PARCIAL – faltam `DVT0W/DVT1W/DVT2W`, `DVT0/DVT1/DVT2` completos com `factor1/sqrt(Xdep)`, `DVTP0..5`, `LPE0/LPEB` |
| Vth – pocket implant (3931) | **não portado** | NÃO – `DVTP0..5` declarados mas sem equação |
| Vth – poly depletion (3984) | `polyDepletion` definido mas **não chamado** (NGATE não portado) | NÃO |
| Vth – Vgsteff (4010) | Vgsteff via Abulk/Esat | PARCIAL |
| Vth – Voff / NFACTOR / MINV | portados | OK (MINV parcial) |
| n (subthreshold ideality, 3906) | n = 1 + NFACTOR | PARCIAL – faltam `CDSC/CDSCB/CDSCD/CIT` |
| Effective channel geometry (4079) | leff/weff | PARCIAL – faltam `XW`, `DWG/DWB`, `LL/LW/...`, `GEOMOD` instância |
| Abulk (4127) | Abulk0 + AGS | PARCIAL – faltam `A1/A2`, `PDIBLC1/2/B` completos, `K3B` termo `LPE0/LPEB` |
| Mobility (4197-4357) | MOBMOD 0..7, high-k | PARCIAL – faltam `UA1/UB1/UC1` (mobilidade 2ª ordem com Vds), `UP/LP`, `EU` (parcial), `UCS/UCSTE`, `RDSW` na mobilidade, `RDSWMIN/RSWMIN/RDWMIN` |
| Vdsat (4357) | Vdsat = EsatL*Vgst2Vtm/(...) | OK |
| Vdseff (4466) | VdseffCV | OK (DC path usa vds direto) |
| Velocity overshoot (4513) | **não portado** | NÃO – `VTL/LC/XN` declarados mas sem equação |
| Vasat (4556) | vsat temp | PARCIAL |
| Idl (4581) | Id = beta*f*lame | OK |
| Degradation pocket (4638) | **não portado** | NÃO |
| VACLM (4651) | PCLM portado (parcial) | PARCIAL – faltam `PDIBLC1/2/B` completos, `PCLM` completo, `PVAG` portado |
| VADIBL (4705) | DIBL simplificado | PARCIAL – faltam `DSUB/DROUT/FPROUT` completos (portados parcialmente), `PSCBE1/2` portado, `ETA0/ETAB` |
| VADITS (4764) | **não portado** | NÃO – `PDITS/PDITSL/PDITSD` declarados |
| VASCBE (4791) | portado | OK |
| Add DIBL/DITS/CLM to Ids (4815-4832) | DIBL gds factor, CLM lambda | PARCIAL |
| Substrate current (4848) | ALPHA0/ALPHA1 | PARCIAL – faltam `BETA0` completo, `AGS` na modulação |
| GIDL (B4p82) | AGIDL/BGIDL/CGIDL/EGIDL/FGIDL/KGIDL/RGIDL | PARCIAL – faltam `RGIDL` (resistência série no GIDL), `GIDLCLAMP` (portado) |
| GISL | AGISL/... | PARCIAL – mesmo acima |
| Gate tunneling IGC/IGB (4.8.2) | AIGC/BIGC/CIGC/AIGS/... | PARCIAL – faltam `PIGCD`, `POXEDGE`, `NTOX`, `EIGBINV`, `NIGBACC/NIGBINV`, `BIGBACC/BIGBINV/CIGBACC/CIGBINV`, `AECHVB/BECHVB` |
| Poly depletion (3984) | função definida, não chamada | NÃO (integração pendente) |
| Velocity overshoot (4513) | **não portado** | NÃO |

---

## 3. Geometria (NumFingerDiff / PAeffGeo / RdseffGeo)

| Item | Status |
|---|---|
| `NumFingerDiff` (min 0 e 1) | OK |
| `PAeffGeo` (casos 0..10) | OK |
| `PERMOD` (strip gate-edge do perímetro) | OK |
| `GEOMOD` instância (honrar instância, não só modelo) | NÃO – usa `model_->GEOMOD` |
| `RdseffGeo` | NÃO |
| `RdsEndIso` / `RdsEndSha` | NÃO (aproximação linear em `evaluateBase`) |
| `DMCGT` (temperatura) | NÃO |
| `DWJ` (junction width offset) | NÃO |
| `XGW` (gate width) no `RGATEMOD` | OK (caso 1) |
| `RGEOMOD` | NÃO |

---

## 4. Temperatura (`updateTemperature4p82_` vs. `scaleTo`)

| Item Xyce | e-XSpice | Status |
|---|---|---|
| UTE (mobilidade) | t.u0 = U0*pow(TRatio,UTE) | OK |
| VSAT (AT) | t.vsat *= 1+AT*(TRatio-1) | OK (no device) |
| Junction sat currents (JSS/JSD/JSWS/...) | t.jctS/D | OK |
| Junction caps (CJS/CJD/...) | t.jctS/D.cap* | OK |
| KT1/KT2 | portados | OK |
| KT1L | portado | OK |
| K3/K3B temp | não | PARCIAL |
| AT (vsat temp) | aplicado | OK |
| PRT (Rds temp) | aplicado | OK |
| UCSTE | não | NÃO |
| TVFBSDOFF / TVOFF / TNFACTOR / TETA0 | não | NÃO |
| TBGASUB/TBGBSUB (bandgap) | usa egfet simplificado | PARCIAL |
| EG0 (high-k) | parcial | PARCIAL |
| NDEP/NSUB temp | não | PARCIAL |
| Self-heating (TEMPMOD) | equilíbrio estático T=Tamb+RTH0*P | PARCIAL – dinâmico (CTH0) não portado |
| SHMOD/RTH0/CTH0/WTH0 | RTH0 usado; CTH0/WTH0/SHMOD não | PARCIAL |

---

## 5. Mobilidade

| Item | Status |
|---|---|
| MOBMOD 0..7 | OK |
| UA/UB/UC (temp escalados) | OK |
| UA1/UB1/UC1 (2ª ordem com Vds) | NÃO |
| UP/LP (mobilidade canal curto) | NÃO |
| EU (expoente) | PARCIAL |
| UCS/UCSTE | PARCIAL |
| UTE | OK |
| LOD (KU0/KVSAT/KVTH0/...) | OK |
| WPE (KU0WE/KVTH0WE/K2WE) | OK |
| Stress (SA/SB/SD/SCA/SCB/SCC) | OK |
| High-k (MTRLMOD/MTRLCOMPATMOD) | PARCIAL |

---

## 6. Vth / body effect / DIBL / CLM

| Item | Status |
|---|---|
| K1/K2 (deriveBodyEffectCoeffs) | OK |
| K3/K3B | OK |
| KETA | OK |
| DVT0/DVT1/DVT2 (simplificado) | PARCIAL |
| DVT0W/DVT1W/DVT2W | NÃO |
| DVTP0..5 | NÃO |
| LPE0/LPEB | NÃO |
| VTH0 L/W dep | NÃO |
| CDSC/CDSCB/CDSCD/CIT | NÃO |
| A0/A1/A2 | A0 OK; A1/A2 NÃO |
| AGS | PARCIAL (portado no Abulk, mas não completo) |
| PDITS/PDITSD/PDITSL | NÃO |
| DELTA | NÃO |
| DSUB/DROUT/FPROUT | PARCIAL (portados parcialmente) |
| ETA0/ETAB | OK |
| PDIBLC1/PDIBLC2/PDIBLCB | PARCIAL |
| PCLM | PARCIAL |
| PVAG | PARCIAL |
| PSCBE1/PSCBE2 | OK |
| DIBL (eta*theta*vds) | OK |

---

## 7. Junctions (saturation + C-V)

| Item | Status |
|---|---|
| Saturation bottom/sidewall/gate-edge | OK |
| Temp scaling (JSS/JSD/JSWS/...) | OK |
| M scaling | OK |
| Junction tunneling (JTSS/...) | PARCIAL (simplificado) |
| Breakdown (BVD/BVS) | PARCIAL (simplificado) |
| IJTHDFWD/SFWD/DREV/SREV | NÃO |
| XJBVD/XJBVS | NÃO |
| VTSS/VTSD/... por termo | PARCIAL |
| NJTS exponents por termo | PARCIAL |
| XTSS/... temp exponents | PARCIAL |
| TNJTS/... | PARCIAL |
| C-V: capMod 0 (Meyer) | PARCIAL (aproximação) |
| C-V: capMod 1/2 (ctmCharges) | OK |
| C-V: capMod 3 | NÃO |
| CVCHARGEMOD 0/1/2 | PARCIAL |
| Overlap (DLC/DWC, CGDO/CGSO/CGBO) | OK |
| Overlap smoothing (CKAPPAS/CKAPPAD/CF) | OK |
| Fringe (CLC/CLE) | OK |
| VFB/VFBCV/VOFFCV/VOFFCVL/MINVCV/ACDE | PARCIAL |
| Poly depletion (NGATE) | NÃO (função existe, não integrada) |

---

## 8. Nós internos (RDSMOD/RGATEMOD/RBODYMOD/RGEOMOD)

| Item | Status |
|---|---|
| RDSMOD (bias-dependent dentro do device) | OK |
| RSH (sheet) | OK |
| RGATEMOD == 1 | OK |
| RGATEMOD == 2/3 | WARN (fallback 0) |
| RBODYMOD == 1 | OK |
| RBODYMOD == 2 | WARN (fallback 1) |
| RGEOMOD | NÃO |
| GEOMOD instância | NÃO |

---

## 9. C-V / cargas

| Item | Status |
|---|---|
| Charges (qg/qd/qb) | OK |
| C-V Jacobian (diferenças centrais) | OK |
| XPART (0/0.5/1) | OK |
| Poly depletion | NÃO (integrado) |
| High-k / metal gate (PHIG/EPSRGATE/...) | PARCIAL (params declarados, equações parciais) |
| EOT scaling (EOT/VDDEOT/...) | NÃO |

---

## 10. Noise

| Item | Status |
|---|---|
| Thermal (RNOIA/RNOIB/RNOIC) | OK |
| Gate/source thermal | OK |
| Flicker (KF/AF/EF) | OK |
| Flicker (NOIA/NOIB/NOIC) | OK |
| FNOIMOD/TNOIMOD | PARCIAL (FNOIMOD 0/1; TNOIMOD não) |
| eval1ovFNoise | PARCIAL (simplificado) |
| getNumNoiseSources | PARCIAL (fixo 2) |
| TNOIA/TNOIB/TNOIC/NTNOI/LINTNOI | NÃO |
| Thermal noise completo (canal + gate + body) | PARCIAL |

---

## 11. Gate/substrate tunneling

| Item | Status |
|---|---|
| AIGC/BIGC/CIGC (gate current) | PARCIAL (simplificado) |
| AIGS/AIGD | PARCIAL |
| AIGBACC/AIGBINV | PARCIAL |
| PIGCD/POXEDGE/NTOX/EIGBINV | NÃO |
| NIGBACC/NIGBINV | NÃO |
| BIGBACC/BIGBINV/CIGBACC/CIGBINV | NÃO |
| IDOVVDS | NÃO |

---

## 12. Body diode / breakdown

| Item | Status |
|---|---|
| Diode básico (saturation) | OK |
| Breakdown (BVD/BVS) | PARCIAL |
| IJTH*/XJB* | NÃO |
| BV/IBV | NÃO |
| Junction temp exponents por termo | PARCIAL |

---

## 13. Binagem / corners

| Item | Status |
|---|---|
| LMIN/LMAX/WMIN/WMAX | OK (honrados) |
| Binagem de parâmetros (L/W/NF) | NÃO (Model::processParams4p82_ do Xyce faz binagem; e-XSpice só checa range) |

---

## 14. High-k / metal gate / poly depletion

| Item | Status |
|---|---|
| PHIG/EPSRGATE/EASUB/EPSRSUB/NI0SUB/BG0SUB/TBGASUB/TBGBSUB | PARCIAL |
| NDEP/NSUB/NSD | PARCIAL |
| POLYDEPLETION | NÃO (função existe, não integrada) |
| EOT/VDDEOT/TEMPEOT/LEFFEOT/WEFFEOT/ADOS/BDOS | NÃO |

---

## 15. NQS (TRNQSMOD/ACNQSMOD)

| Item | Status |
|---|---|
| TRNQSMOD | WARN (desabilitado) |
| ACNQSMOD | WARN (desabilitado) |
| TAUN/ELM/XRCRG1/XRCRG2/QDEF | NÃO (QDEF declarado, não usado) |

Motivo: requer infraestrutura de estado DAE (state vectors, Master::updateState) não disponível no e-XSpice.

---

## 16. Código latente / morto encontrado na verificação

- `MosfetBSIM4Device::updatePrimaryState()`, `auxChargeCalculations()`, `setupCapacitors()` estavam **declarados mas não definidos nem chamados**. Removidos do header para evitar símbolos indefinidos latentes.
- `std::fprintf("[dbg] ...")` em `evaluateBase` imprimia em toda avaliação (não guardado por `debug_`). Removido.
- `eval1ovFNoise` e `getNumNoiseSources` (definidos na P2) são funcionais mas simplificados.
- `polyDepletion` definido mas não chamado no fluxo principal (NGATE não portado).
- Parâmetros verdadeiramente não referenciados em nenhum arquivo: `AECHVB`, `BECHVB`, `cgd`, `DTEMP`, `EM`, `ketac`, `RS0`, `thetaRout`, `VFBSD`, `VGSTEFFVTH`, `wlc`, `wlr`. A maioria é legado/alias; manter por compatibilidade de netlist.

---

## 17. Resumo quantitativo (atualizado)

- Funções Xyce Instance relevantes: ~30. e-XSpice implementa ~22 (73%).
- Blocos físicos de `updateIntermediateVars4p82_`: ~25. Portados: ~16 (64%), parciais: ~6, não portados: ~3 (pocket implant, velocity overshoot, VADITS) + integração pendente (poly depletion, EOT).
- Parâmetros declarados: ~485. Usados em physics: ~126 (Device.cpp) + ~359 (Model.cpp/Mosfet.cpp). Não referenciados: 12.
- Modos desabilitados com warning: RGATEMOD 2/3, RBODYMOD 2, TRNQSMOD, ACNQSMOD, TEMPMOD dinâmico, NQS.

---

## 18. Prioridade sugerida (atualizada)

1. **Alta (física core):** DVT*/CDSC/CIT completo, CLM completo (A1/A2, PDIBLC*, PCLM, PVAG), poly depletion integrado, mobilidade UA1/UB1/UC1 + UP/LP, velocity overshoot, gate tunneling completo (PIGCD/POXEDGE), junction IJTH*/XJB*.
2. **Média (geometria):** RdseffGeo/RdsEndIso/Sha, GEOMOD instância, RGEOMOD, DMCGT, DWJ.
3. **Média (temperatura):** UCSTE, TVFBSDOFF/TVOFF/TNFACTOR/TETA0, TBGASUB/TBGBSUB, bandgap completo, CTH0 dinâmico.
4. **Baixa (avançado):** NQS (requer estado DAE), binagem de parâmetros, EOT scaling, high-k completo, capMod 3, CVCHARGEMOD completo.
5. **Infra:** definir `eval1ovFNoise` completo (já parcial), `getNumNoiseSources` dinâmico, remover código morto.

---

*Gerado por verificação contra `external_codes/Xyce/OpenModels/N_DEV_MOSFET_B4{,p82}.C`.*
