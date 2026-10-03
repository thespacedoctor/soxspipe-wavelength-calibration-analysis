#!/usr/bin/env python3
"""
Generate a telluric absorption line list from HITRAN.

Default range: 8000-8400 A  (a dense H2O band; O2 contributes essentially nothing here)

Requires:  pip install hitran-api numpy
Needs outbound access to hitran.org (the HITRAN API is rate-limited: a few
queries per day for anonymous users, so the local cache below matters).

Usage:
    python telluric_lines_from_hitran.py
    python telluric_lines_from_hitran.py --wmin 7800 --wmax 8400 --min-intensity 1e-25
"""

import argparse
import os

import numpy as np
from hapi import db_begin, fetch, getColumn, tableList

# ---------------------------------------------------------------- conversions

def vac_to_air(w):
    """Vacuum -> air wavelength in Angstroms (Ciddor 1996; the VALD/NIST standard).

    Validated against Ca II 8542.09 A (0.003 A) and H-alpha 6562.82 A (0.02 A).
    """
    s2 = (1e4 / np.asarray(w, dtype=float)) ** 2
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2)
    return w / n


def wavenumber_to_wavelength(nu):
    """HITRAN wavenumber (cm^-1, vacuum) -> vacuum wavelength (Angstroms)."""
    return 1e8 / np.asarray(nu, dtype=float)


# ---------------------------------------------------------------- main routine

MOLECULES = {
    # name: (HITRAN molecule id, isotopologue id, global iso id)
    "H2O": (1, 1, 1),
    "O2": (7, 1, 36),
}


def build_line_list(wmin_a, wmax_a, min_intensity, air=True, cache="hitran_cache"):
    # Angstroms (air-ish input) -> wavenumbers. Pad the window by 1 A so that
    # lines near the edge are not clipped by the vac/air shift (~2 A at 8000 A).
    nu_min = 1e8 / (wmax_a + 5.0)
    nu_max = 1e8 / (wmin_a - 5.0)

    os.makedirs(cache, exist_ok=True)
    db_begin(cache)

    rows = []
    for name, (mol_id, iso_id, _) in MOLECULES.items():
        table = f"{name}_{int(wmin_a)}_{int(wmax_a)}"
        if table not in tableList():
            print(f"fetching {name} from HITRAN ({nu_min:.1f}-{nu_max:.1f} cm-1) ...")
            fetch(table, mol_id, iso_id, nu_min, nu_max)
        else:
            print(f"using cached table {table}")

        nu = np.array(getColumn(table, "nu"))
        sw = np.array(getColumn(table, "sw"))  # line intensity, cm-1/(molec cm-2)

        keep = sw >= min_intensity
        for n, s in zip(nu[keep], sw[keep]):
            rows.append((n, s, name))

    if not rows:
        raise SystemExit("no lines returned -- check the range or lower --min-intensity")

    rows.sort()
    nu = np.array([r[0] for r in rows])
    sw = np.array([r[1] for r in rows])
    species = [r[2] for r in rows]

    wl = wavenumber_to_wavelength(nu)
    if air:
        wl = vac_to_air(wl)

    # trim back to the requested window after the shift
    m = (wl >= wmin_a) & (wl <= wmax_a)
    order = np.argsort(wl[m])
    return wl[m][order], sw[m][order], [s for s, k in zip(species, m) if k]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--wmin", type=float, default=8000.0, help="min wavelength, Angstroms")
    p.add_argument("--wmax", type=float, default=8400.0, help="max wavelength, Angstroms")
    p.add_argument("--min-intensity", type=float, default=1e-25,
                   help="HITRAN sw cutoff. 1e-25 gives the strong lines; "
                        "1e-26 ~ several hundred; 1e-27 ~ thousands.")
    p.add_argument("--vacuum", action="store_true",
                   help="output vacuum wavelengths (default is air, matching "
                        "the Catania atlas and most optical work)")
    p.add_argument("--out", default="telluric_8000_8400.py")
    args = p.parse_args()

    wl, sw, species = build_line_list(
        args.wmin, args.wmax, args.min_intensity, air=not args.vacuum
    )

    frame = "vacuum" if args.vacuum else "air"
    print(f"\n{len(wl)} lines, {args.wmin:.0f}-{args.wmax:.0f} A ({frame})")

    with open(args.out, "w") as f:
        f.write(f"# telluric absorption lines, {args.wmin:.0f}-{args.wmax:.0f} A ({frame})\n")
        f.write(f"# HITRAN, intensity cutoff sw >= {args.min_intensity:g}\n")
        f.write("telluric_lines = [\n    ")
        f.write(", ".join(f"{w:.2f}" for w in wl))
        f.write("\n]\n")
    print(f"wrote {args.out}")

    # also dump wavelength / intensity / species for inspection
    with open(args.out.replace(".py", ".txt"), "w") as f:
        f.write("# wavelength_A  intensity  species\n")
        for w, s, sp in zip(wl, sw, species):
            f.write(f"{w:10.3f}  {s:.3e}  {sp}\n")


if __name__ == "__main__":
    main()
