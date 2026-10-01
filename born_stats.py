import sys
import numpy as np
from itertools import islice

BAR_TO_GPA = 1e-4
ROOT   = '/scratch/hpcl-mat269/lls/PRODUCTION_DATA'
DELTA  = 1.0e-3
VFIELD = 'c_vor[1]'
VORODUMP = 'als_{T}.dump'
STRIDE, MAX_FRAMES = 10, 1000

def _iter_blocks(path, keep):
    with open(path) as f:
        while True:
            h = [f.readline() for _ in range(9)]
            if not h[0]:
                break
            n = int(h[3])
            arr = np.array([f.readline().split() for _ in range(n)], dtype=float)
            yield arr[:, keep]

def load_omega(als, vol_field=VFIELD, stride=STRIDE, max_frames=MAX_FRAMES):
    with open(als) as f:
        lines = f.readlines()
    hdr = next(l for l in lines if l.startswith('ITEM: ATOMS')).split()[2:]
    ci, ct, cv = hdr.index('id'), hdr.index('type'), hdr.index(vol_field)
    acc = None; cnt = 0; types = None; i = 0; frame_idx = 0
    while i < len(lines):
        n = int(lines[i+3])
        if frame_idx % stride == 0:
            arr = np.array([ln.split() for ln in lines[i+9:i+9+n]], dtype=float)
            arr = arr[np.argsort(arr[:, ci])]
            if acc is None:
                acc = np.zeros(n); types = arr[:, ct].astype(int)
            acc += arr[:, cv]; cnt += 1
            if cnt >= max_frames:
                break
        i += 9 + n; frame_idx += 1
    return acc/cnt, types

def elastic_species(bornpa, omega, types, delta, stride=STRIDE, max_frames=MAX_FRAMES):
    """
    Per-atom affine (Born) elastic constants from axial strains, species-averaged
        C11 = <diag>_i , C12 = <off-diag>_i
        B   = (C11 + 2 C12)/3           bulk modulus
        G   = (C11 - C12)/2             shear modulus, affine isotropic
        nu  = C12/(C11 + C12)           Poisson ratio = (3B-2G)/(2(3B+G))
        E   = 2 G (1 + nu)              Young's modulus
    """
    gen = _iter_blocks(bornpa, [2, 3, 4])
    cu, zr = types == 2, types == 1
    rows = {k: [] for k in ("C11_Cu","C12_Cu","B_Cu","G_Cu","nu_Cu","E_Cu",
                            "C11_Zr","C12_Zr","B_Zr","G_Zr","nu_Zr","E_Zr",
                            "B_all","G_all","nu_all","E_all")}
    frame_idx = 0
    while True:
        b = list(islice(gen, 6))
        if len(b) < 6:
            break
        if frame_idx % stride == 0:
            D = np.stack([(b[0]-b[1])/(2*delta),          # Omega*dSxx,yy,zz / d eps_x
                          (b[2]-b[3])/(2*delta),          #                    / d eps_y
                          (b[4]-b[5])/(2*delta)])         #                    / d eps_z
            diag = (D[0, :, 0] + D[1, :, 1] + D[2, :, 2]) / 3.0
            off  = (D[0, :, 1] + D[1, :, 0] + D[0, :, 2] +
                    D[2, :, 0] + D[1, :, 2] + D[2, :, 1]) / 6.0
            C11_i = diag / omega * BAR_TO_GPA
            C12_i = off  / omega * BAR_TO_GPA
            for lab, m in (("Cu", cu), ("Zr", zr), ("all", slice(None))):
                c11, c12 = C11_i[m].mean(), C12_i[m].mean()
                B = (c11 + 2*c12)/3.0
                G = (c11 - c12)/2.0
                nu = c12/(c11 + c12)
                E = 2*G*(1 + nu)
                if lab != "all":
                    rows[f"C11_{lab}"].append(c11); rows[f"C12_{lab}"].append(c12)
                rows[f"B_{lab}"].append(B); rows[f"G_{lab}"].append(G)
                rows[f"nu_{lab}"].append(nu); rows[f"E_{lab}"].append(E)
            if len(rows["B_Cu"]) >= max_frames:
                break
        frame_idx += 1
    return {k: np.array(v) for k, v in rows.items()}

def run(T):
    omega, types = load_omega(f'{ROOT}/{VORODUMP.format(T=T)}')
    r = elastic_species(f'{ROOT}/bornpa_{T}.dump', omega, types, DELTA)
    def ms(k):
        a = r[k]; return a.mean(), a.std(ddof=1)/np.sqrt(len(a))
    (K_Cu, dK_Cu), (K_Zr, dK_Zr) = ms("B_Cu"), ms("B_Zr")
    (G_Cu, dG_Cu), (G_Zr, dG_Zr) = ms("G_Cu"), ms("G_Zr")
    (nu_Cu, dnu_Cu), (nu_Zr, dnu_Zr) = ms("nu_Cu"), ms("nu_Zr")
    (B_all, dB_all), (G_all, dG_all) = ms("B_all"), ms("G_all")
    (nu_all, dnu_all), (E_all, dE_all) = ms("nu_all"), ms("E_all")
    np.savez(f'{ROOT}/Bs_{T}.npz',
             K_Cu=K_Cu, K_Zr=K_Zr, dK_Cu=dK_Cu, dK_Zr=dK_Zr,
             K_Cu_frames=r["B_Cu"], K_Zr_frames=r["B_Zr"],
             G_Cu=G_Cu, G_Zr=G_Zr, dG_Cu=dG_Cu, dG_Zr=dG_Zr,
             G_Cu_frames=r["G_Cu"], G_Zr_frames=r["G_Zr"],
             nu_Cu=nu_Cu, nu_Zr=nu_Zr, dnu_Cu=dnu_Cu, dnu_Zr=dnu_Zr,
             C11_Cu_frames=r["C11_Cu"], C12_Cu_frames=r["C12_Cu"],
             C11_Zr_frames=r["C11_Zr"], C12_Zr_frames=r["C12_Zr"],
             B_all=B_all, dB_all=dB_all, G_all=G_all, dG_all=dG_all,
             nu_all=nu_all, dnu_all=dnu_all, E_all=E_all, dE_all=dE_all,
             B_all_frames=r["B_all"], G_all_frames=r["G_all"],
             nu_all_frames=r["nu_all"], E_all_frames=r["E_all"], T=T)
    print(f'{T}  B_Cu={K_Cu:6.2f}±{dK_Cu:.2f}  B_Zr={K_Zr:6.2f}±{dK_Zr:.2f} | '
          f'G_all={G_all:6.2f}±{dG_all:.2f}  B_all={B_all:6.2f}  '
          f'nu_all={nu_all:.3f}±{dnu_all:.3f}  E_all={E_all:6.2f} GPa '
          f'(N={len(r["B_Cu"])})', flush=True)

if __name__ == '__main__':
    run(int(sys.argv[1]))
