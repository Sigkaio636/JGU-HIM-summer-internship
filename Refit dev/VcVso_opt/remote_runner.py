import yaml
from pathlib import Path
import subprocess as subp
from dataclasses import dataclass

import re
import csv

# <> - <> - <> - <> - <>
#   AUXILIARY METHODS
# <> - <> - <> - <> - <>

EBG = 2.467000
EBE = 1.990390


@dataclass
class DataRow:
    iter: int
    Vc: float
    Vso: float
    Eg: float
    Ee: float
    L: float
    Dx: float
    lr: float


def create_inpfile_VcVso(vc, vso, nameidx):
    template = """\
    197.32705359 1.43996518 20.736    * csts: hc,e2,hm
    4.0 3.0 2 1 0 1         * Masses, charges, spins: Ac,Af (in amu),Zc,Zf,iSc,iSf
    1 1                     * c-f potential: ntypo,Njpo
    {Vc:.8f} 2.39 0.68 {Vso:.8f} 2.39 0.68 2.39        *j,l,Vp,rp,ap,VLS,rls,als,rC
    1                       * nb of bound states
    1 {s:1d} 1           * bound states: NR,J,l
    20000 100 1e-6          * radial mesh:Nru,rmu,eps
    1500 3.01 0.01          * Energies: NE,E0,hE
    0                       * lec 0:E, 10:wf 1st state, 100:phaseshift, 110 bin
    Li7comp
    """

    out_dir = Path("inputs")
    out_dir.mkdir(exist_ok=True)

    state = 3 if nameidx[5] == "G" else 1
    content = template.format(Vc=vc, Vso=vso, s=state)
    filename = out_dir / f"{nameidx}_input.dat"
    filename.write_text(content)


def execute_Boscos(nameidx):
    with (
        open(f"inputs/{nameidx}_input.dat", "r") as infile,
        open(f"outputs/{nameidx}_output.txt", "w") as outfile,
    ):
        result = subp.run(
            ["../Boscos.out"],
            stdin=infile,
            stdout=outfile,
            stderr=subp.PIPE,
            check=True,
        )
        # print("STDERR:", result.stderr)
        # print("STDOUT:", result.stdout)


def read_outputed_E(nameidx):
    with open(f"outputs/{nameidx}_output.txt", "r") as f:
        content = f.read()

        match = re.search(r"E=\s*([-+]?\d+\.\d+)\s+MeV", content)
        if match:
            energy = float(match.group(1))
        else:
            energy = None

    assert not energy is None, f"Not bounded state found at iter {nameidx}."
    return energy


def pipeline_Boscos(vc, vso, nameidx):
    create_inpfile_VcVso(vc, vso, nameidx)
    execute_Boscos(nameidx)
    return read_outputed_E(nameidx)


def evaluate_L_and_grad_and_Ebs(name_num, vc, vso, h=1e-4):
    nameidx = f"{name_num:04d}"
    # offset:
    # 0:center, 1:+hx, 2:-hx, 3:+hy, 4:-hy
    Bg = [0.0] * 5
    Be = [0.0] * 5

    Bg[0] = pipeline_Boscos(vc, vso, nameidx + "_G_0")
    Bg[1] = pipeline_Boscos(vc + h, vso, nameidx + "_G_1")
    Bg[2] = pipeline_Boscos(vc - h, vso, nameidx + "_G_2")
    Bg[3] = pipeline_Boscos(vc, vso + h, nameidx + "_G_3")
    Bg[4] = pipeline_Boscos(vc, vso - h, nameidx + "_G_4")

    Be[0] = pipeline_Boscos(vc, vso, nameidx + "_E_0")
    Be[1] = pipeline_Boscos(vc + h, vso, nameidx + "_E_1")
    Be[2] = pipeline_Boscos(vc - h, vso, nameidx + "_E_2")
    Be[3] = pipeline_Boscos(vc, vso + h, nameidx + "_E_3")
    Be[4] = pipeline_Boscos(vc, vso - h, nameidx + "_E_4")

    Lossfunc = 0.25 * ((Bg[0] + EBG) ** 2 + (Be[0] + EBE) ** 2)
    # print(Bg, Be)

    Bg_vc = (Bg[1] - Bg[2]) / (2 * h)
    Bg_vso = (Bg[3] - Bg[4]) / (2 * h)
    Be_vc = (Be[1] - Be[2]) / (2 * h)
    Be_vso = (Be[3] - Be[4]) / (2 * h)
    Loss_vc = 0.5 * (Bg[0] + EBG) * Bg_vc + 0.5 * (Be[0] + EBE) * Be_vc
    Loss_vso = 0.5 * (Bg[0] + EBG) * Bg_vso + 0.5 * (Be[0] + EBE) * Be_vso

    return Lossfunc, Loss_vc, Loss_vso, Bg[0], Be[0]


def append_to_csv(row):
    with open("run_log.csv", "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [row.iter, row.Vc, row.Vso, row.Eg, row.Ee, row.L, row.Dx, row.lr]
        )


# <> - <> - <> - <> - <>
#   MAIN METHODS
# <> - <> - <> - <> - <>


def main():

    # Initialize first point
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)
    vp0 = config.get("sweep", {}).get("x0", [-80, -10])

    L, L_vc, L_vso, EbG, EbE = evaluate_L_and_grad_and_Ebs(0, vp0[0], vp0[1])
    lr = 12.73086336173236 / 1
    append_to_csv(DataRow(0, vp0[0], vp0[1], EbG, EbE, L, 0, 0))

    print(L, L_vc, L_vso, lr)

    # Iterate the gradient descend
    vp_l = vp0
    idx = 1
    while abs(L_vc) + abs(L_vso) > 1e-6:
        vp_n = (vp_l[0] - lr * L_vc, vp_l[1] - lr * L_vso)
        L, L_vc, L_vso, EbG, EbE = evaluate_L_and_grad_and_Ebs(idx, vp_n[0], vp_n[1])
        Dx = abs(vp_n[0] - vp_l[0]) + abs(vp_n[1] - vp_l[1])
        append_to_csv(DataRow(idx, vp_n[0], vp_n[1], EbG, EbE, L, Dx, lr))

        idx += 1
        vp_l = vp_n


main()
