# AutoDock4 Parameter Profiles

ADSP organizes AutoDock4 force-field parameter profiles into isolated, auditable profiles with cryptographic verification.

## Available Profiles

- **`ad4_standard_4.2`**: Standard non-specialized AutoDock4.2 parameter set (`AD4_parameters.dat`).
  - Supported atom types: 35 standard types including H, HD, C, A, N, NA, OA, OS, P, SA, S, Cl, Br, I, F, Mg, Ca, Mn, Fe, Zn.
  - SHA-256: `625DE5779B914382E21A135C776EFBC02B4221085BD0280118D103CCDD93EA7C`
- **`ad4_1_bound`**: AutoDock4.1 bound-state parameter set (`AD4.1_bound.dat`).
  - Reference: Huey et al. (2007), J Comput Chem 28:1145-1152.
  - SHA-256: `6B98F7AB508F4882801938F8CED1C0BF38096496155A8005BAF941A201781CE8`

## Profile Isolation

Profiles are never merged or blended implicitly. The user or configuration explicitly selects the parameter set for GPF and AutoGrid4 generation.
