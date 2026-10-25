# BSIM4 Wrapper Architecture: Berkeley C Code Integration

## Status

The Berkeley BSIM4 4.8.3 C code is vendored under `src/models/bsim4/code/` (29 files,
~100K lines). It is written against the SPICE device interface (`CKTcircuit`,
`GENinstance`, `GENmodel`, `SMPmatrix`, `SMPaddElt`, etc.) and requires the
SPICE framework headers (`spice.h`, `gendefs.h`, `devdefs.h`, `suffix.h`,
`cktdefs.h`, `ifsim.h`, `sperror.h`, `util.h`).

A minimal shim of these headers exists in `src/models/external/SpiceContext.h`,
but the Berkeley code uses many additional SPICE internals (state vectors,
matrix pointers, node lists, mode flags, etc.) that are not yet shimmed.

## What is needed for full integration

### 1. Complete SPICE compatibility layer

The Berkeley BSIM4 code accesses the following SPICE internals:

- `CKTcircuit`: CKTtemp, CKTgmin, CKTmode, CKTrhs, CKTstate0, CKTstate1,
  CKTmatrix, CKTnodes, CKTnumNodes, CKTnode, CKTcurPtr, CKTorder, CKTstep,
  CKTdelta, CKTdeltaOld, CKTdeltaOld2, CKTdeltaOld3, CKTdeltaOld4,
  CKTdeltaOld5, CKTdeltaOld6, CKTdeltaOld7, CKTdeltaOld8, CKTdeltaOld9,
  CKTdeltaOld10, CKTdeltaOld11, CKTdeltaOld12, CKTdeltaOld13,
  CKTdeltaOld14, CKTdeltaOld15, CKTdeltaOld16, CKTdeltaOld17,
  CKTdeltaOld18, CKTdeltaOld19, CKTdeltaOld20, ...
- `GENinstance`: GENnumNodes, GENnode, GENstate, GENparam, GENmodel
- `GENmodel`: GENmodName, GENmodType, GENinstances, GENnextModel
- `SMPmatrix`: SMPsize, SMPmatrix, SMPmatrixSize, SMPmatrixPtr
- `SMPaddElt`: returns a pointer to the matrix element
- `CKTnode`: CKTnodeNumber, CKTnodeName
- `IFparm`, `IFvalue`, `IFuid`: parameter tables
- `GENERIC` macros, `DEV_DEFAULT`, etc.

This is essentially the entire SPICE device interface. Implementing it
requires building a complete SPICE compatibility layer.

### 2. The adapter

The adapter (`Bsim4Adapter : SpiceDeviceAdapter`) would:
1. In `registerJacobianPattern`: run `BSIM4setup` to discover the (row, col)
   pairs the device will stamp, and declare them in the e-XSpice pattern.
2. In `loadF`/`loadQ`/`loadJacobianF`/`loadJacobianQ`/`loadB`: build the
   shimmed `CKTcircuit`, call the Berkeley load function (`BSIM4load`), and
   let the `SMPaddElt` shim translate the stamps into the pre-resolved slots.
3. In `loadNoiseSources`: call `BSIM4noise`.
4. In temperature updates: call `BSIM4temp`.

### 3. The wiring

- Compile the Berkeley C code as a static library (with the SPICE shim).
- Link it into `e-XSpice_lib`.
- Register the `Bsim4Adapter` in the device factory.

## Effort estimate

This is a multi-month project. The SPICE compatibility layer alone is
~10K lines of code. The Berkeley BSIM4 code is ~100K lines and uses many
SPICE internals.

## Recommended approach

1. **Start with the SPICE shim**: implement `CKTcircuit`, `SMPmatrix`,
   `SMPaddElt`, and the node list. This is the foundation.
2. **Implement the adapter**: create `Bsim4Adapter` that wraps the Berkeley
   code. Start with DC evaluation (BSIM4load) and add AC/noise/charge as
   needed.
3. **Test incrementally**: start with a simple MOSFET circuit and compare
   against the existing BSIM4 implementation.
4. **Add EKV**: follow the same pattern with the EKV C code.

## Current status

- The architecture is in place (`SpiceContext.h`, `SpiceDeviceAdapter.h/cpp`).
- The "level 2" naming has been fixed.
- The Berkeley C code is vendored under `src/models/bsim4/code/`.
- The SPICE compatibility layer is minimal (only the fields the existing
  code needs are implemented).
- The full integration of the Berkeley C code requires building the
  complete SPICE compatibility layer, which is a substantial undertaking.

## Why this is the right architecture

The SPICE device adapter pattern is the same one used by ngspice (which IS
the Berkeley code, integrated) and by Xyce (for its ADMS-generated devices).
The key insight is that the Berkeley code is written against a well-defined
interface (`CKTcircuit`, `GENinstance`, `SMPaddElt`), and the adapter only
needs to provide that interface — the physics is never re-implemented.

The e-XSpice DAE contract is the canonical interface, and the SPICE shim
is an ADAPTER, not the core. This means:
- The analysis engines (`.OP`, `.TRAN`, `.AC`, `.NOISE`, `.TF`, `.PZ`,
  `.SENS`, `.HB`) are device-agnostic: they consume the DAE contract.
- Adding a new model never touches the analysis engines.
- The shim grows by adding fields, never by rewriting the adapter.
- The physics is the Berkeley code's, not a re-implementation.