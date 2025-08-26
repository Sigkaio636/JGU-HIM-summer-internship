import pathlib
import yaml

# Lee el yaml
script_dir = pathlib.Path(__file__).parent.resolve()
with open(f"{script_dir}/config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Extrae Vc_values
Vc_param = config.get("sweep", {}).get("Vc", [-80, 2, -70])  # default por si no está
Vc_values = []
while Vc_param[0] <= Vc_param[2]:
    Vc_values.append(round(Vc_param[0], 10))
    Vc_param[0] += Vc_param[1]

template = """\
197.32705359 1.43996518 20.900796404    * csts: hc,e2,hm
4.002603254 3.016049281 2 1 0 1         * Masses, charges, spins: Ac,Af (in amu),Zc,Zf,iSc,iSf
1 1                     * c-f potential: ntypo,Njpo
{Vc:.8f} 2.39 0.68 0.00 2.39 0.68 2.35        *j,l,Vp,rp,ap,VLS,rls,als,rC
1                       * nb of bound states
1 1 1           * bound states: NR,J,l
20000 100 1e-6          * radial mesh:Nru,rmu,eps
1500 3.01 0.01          * Energies: NE,E0,hE
0                       * lec 0:E, 10:wf 1st state, 100:phaseshift, 110 bin
Li7comp
"""

out_dir = pathlib.Path("inputs")
out_dir.mkdir(exist_ok=True)

for i, vc in enumerate(Vc_values, start=1):
    content = template.format(Vc=vc)
    filename = out_dir / f"{i:03d}_input.dat"
    filename.write_text(content)
    print(f"  → Created {filename}")

"""
Conventions
Eki:

197.32705359 1.43996518 20.900796404    * csts: hc,e2,hm
4.002603254 3.016049281 2 1 0 1         * Masses, charges, spins: Ac,Af (in amu),Zc,Zf,iSc,iSf
1 1                     * c-f potential: ntypo,Njpo
{Vc:.8f} 2.39 0.68 0.00 2.39 0.68 2.33951        *j,l,Vp,rp,ap,VLS,rls,als,rC
1                       * nb of bound states
1 3 1           * bound states: NR,J,l
20000 100 1e-6          * radial mesh:Nru,rmu,eps
1500 3.01 0.01          * Energies: NE,E0,hE
0                       * lec 0:E, 10:wf 1st state, 100:phaseshift, 110 bin
Li7comp

Pierre:

197.32705359 1.43996518 20.736    * csts: hc,e2,hm
4.0 3.0 2 1 0 1         * Masses, charges, spins: Ac,Af (in amu),Zc,Zf,iSc,iSf
1 1                     * c-f potential: ntypo,Njpo
{Vc:.8f} 2.39 0.68 0.00 2.39 0.68 2.39        *j,l,Vp,rp,ap,VLS,rls,als,rC
1                       * nb of bound states
1 3 1           * bound states: NR,J,l
20000 100 1e-6          * radial mesh:Nru,rmu,eps
1500 3.01 0.01          * Energies: NE,E0,hE
0                       * lec 0:E, 10:wf 1st state, 100:phaseshift, 110 bin
Li7comp
"""
