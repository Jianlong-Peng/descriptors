/*=============================================================================
#     FileName: names.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-12-27 12:34:22
#   LastChange: 2014-02-13 20:05:19
#      History:
=============================================================================*/
#include <vector>
#include <string>
#include <openbabel/babelconfig.h>
#include "names.h"
#include "config.h"

using std::vector;
using std::string;

static string str_constitute[] = {"number of atoms",
    "number of C atoms", "relative number of C atoms",
    "number of H atoms", "relative number of H atoms",
    "number of O atoms", "relative number of O atoms",
    "number of N atoms", "relative number of N atoms",
    "number of S atoms", "relative number of S atoms",
    "number of F atoms", "relative number of F atoms",
    "number of Cl atoms", "relative number of Cl atoms",
    "number of Br atoms", "relative number of Br atoms",
    "number of I atoms", "relative number of I atoms",
    "number of P atoms", "relative number of P atoms",
    "molecular weight", "relative molecular weight",
    "number of bonds", "number of single bonds", "relative number of single bonds",
    "number of double bonds", "relative number of double bonds",
    "number of triple bonds", "relative number of triple bonds",
    "number of aromatic bonds", "relative number of aromatic bonds",
    "number of rings", "relative number of rings",
    "number of benzene rings", "relative number of benzene rings",
    "number of H-acceptor sites", "number of H-donor sites"
};
vector<string> DESCRIPTOR_API constitute_names(str_constitute, str_constitute+38);

static string str_electrostatic[] = {
    "max partial charge for H [Zefirov]", "min partial charge for H [Zefirov]",
    "max partial charge for C [Zefirov]", "min partial charge for C [Zefirov]",
    "max partial charge for N [Zefirov]", "min partial charge for N [Zefirov]",
    "max partial charge for O [Zefirov]", "min partial charge for O [Zefirov]",
    "max partial charge for F [Zefirov]", "min partial charge for F [Zefirov]",
    "max partial charge for P [Zefirov]", "min partial charge for P [Zefirov]",
    "max partial charge for S [Zefirov]", "min partial charge for S [Zefirov]",
    "max partial charge for Cl [Zefirov]", "min partial charge for Cl [Zefirov]",
    "max partial charge for Br [Zefirov]", "min partial charge for Br [Zefirov]",
    "max partial charge for I [Zefirov]", "min partial charge for I [Zefirov]",
    "max partial charge (Qmax) [Zefirov]", "min partial charge (Qmin) [Zefirov]",
    "polarity parameter (Qmax-Qmin)", "polarity parameter / square distance",
    "topographic electronic index (all pairs) [Zefirov]", "topographic electronic index (all bonds) [Zefirov]"
};
vector<string> DESCRIPTOR_API electrostatic_names(str_electrostatic,str_electrostatic+26);

static string str_geometrical[] = {
    "moment of inertia A", "moment of inertia B", "moment of inertia C",
    "XY shadow", "ZX shadow", "YZ shadow",
    "XY shadow / XY rectangle", "ZX shadow / ZX rectangle", "YZ shadow / YZ rectangle",
    "molecular volume", "molecular volume / XYZ box",
    "gravitation index (all pairs)", "gravitation index (all bonds)",
    "TPSA"//, "molecular surface area"
};
vector<string> DESCRIPTOR_API geometrical_names(str_geometrical,str_geometrical+14);

static string str_topological[] = {"Wiener index",
    "Randic index (order 0)", "Randic index (order 1)", "Randic index (order 2)", "Randic index (order 3)",
    "Kier&Hall index (order 0)", "Kier&Hall index (order 1)", "Kier&Hall index (order 2)", "Kier&Hall index (order 3)",
    "Kier flexibility index", "Kier shape index (order 1)", "Kier shape index (order 2)", "Kier shape index (order 3)",
    "Balaban index",
    "average information content (order 0)", "information content (order 0)",
    "average structural information content (order 0)", "structural information content (order 0)",
    "average complementary information content (order 0)", "complementary information content (order 0)",
    "average bonding information content (order 0)", "bonding information content (order 0)",
    "average information content (order 1)", "information content (order 1)",
    "average structural information content (order 1)", "structural information content (order 1)",
    "average complementary information content (order 1)", "complementary information content (order 1)",
    "average bonding information content (order 1)", "bonding information content (order 1)",
    "average information content (order 2)", "information content (order 2)",
    "average structural information content (order 2)", "structural information content (order 2)",
    "average complementary information content (order 2)", "complementary information content (order 2)",
    "average bonding information content (order 2)", "bonding information content (order 2)"
};
vector<string> DESCRIPTOR_API topological_names(str_topological,str_topological+38);

static string str_qc[] = {
    "final heat of formation", "final heat of formation / #atoms",
    "No. of occupied electronic levels", "No. of occupied electronic levels / #atoms",
    "HOMO-1 energy", "HOMO energy", "LUMO energy", "LUMO+1 energy",
    "min NA for a C atom", "max NA for a C atom", "aver. NA for a C atom",
    "min NA for a N atom", "max NA for a N atom", "aver. NA for a N atom",
    "min NA for a O atom", "max NA for a O atom", "aver. NA for a O atom",
    "min EA for a C atom", "max EA for a C atom", "aver. EA for a C atom",
    "min EA for a N atom", "max EA for a N atom", "aver. EA for a N atom",
    "min EA for a O atom", "max EA for a O atom", "aver. EA for a O atom",
    "min RA for a C atom", "max RA for a C atom", "aver. RA for a C atom",
    "min RA for a N atom", "max RA for a N atom", "aver. RA for a N atom",
    "min RA for a O atom", "max RA for a O atom", "aver. RA for a O atom",
    "max atomic charge for a C atom [Mulliken]", "min atomic charge for a C atom [Mulliken]",
    "max atomic charge for a N atom [Mulliken]", "min atomic charge for a N atom [Mulliken]",
    "max atomic charge for a O atom [Mulliken]", "min atomic charge for a O atom [Mulliken]",
    "max atomic charge for a H atom [Mulliken]", "min atomic charge for a H atom [Mulliken]",
    "max atomic charge [Mulliken]", "min atomic charge [Mulliken]",
    "Tot. point-charge comp. of the molecular dipole",
    "Tot. hybridization comp. of the molecular dipole",
    "Tot. dipole of the molecule",
    "min atomic orbital electronic population", "max atomic orbital electronic population",
    "max SIGMA-SIGMA bond order", "min SIGMA-SIGMA bond order", 
    "max SIGMA-PI bond order", "min SIGMA-PI bond order", 
    "max PI-PI bond order", "min PI-PI bond order",
    "Tot. of E-N attraction (1-center)", "aver. of E-N attraction (1-center)",
    "Tot. of E-E repulsion (1-center)", "aver. of E-E repulsion (1-center)",
    "Tot. of resonance energy (2-center)", "aver. of resonance energy (2-center)",
    "Tot. of exchange energy (2-center)", "aver. of exchange energy (2-center)",
    "Tot. of E-E repulsion (2-center)", "aver. of E-E repulsion (2-cneter)",
    "Tot. of E-N attraction (2-center)", "aver. of E-N attraction (2-center)",
    "Tot. of N-N repulsion (2-center)", "aver. of N-N repulsion (2-center)",
    "max atomic charge for a C atom [ESP]", "min atomic charge for a C atom [ESP]",
    "max atomic charge for a N atom [ESP]", "min atomic charge for a N atom [ESP]",
    "max atomic charge for a O atom [ESP]", "min atomic charge for a O atom [ESP]",
    "max atomic charge for a H atom [ESP]", "min atomic charge for a H atom [ESP]",
    "max atomic charge [ESP]", "min atomic charge[ESP]"
};
vector<string> DESCRIPTOR_API qc_names(str_qc, str_qc+80);

static string str_physicochemical[] = {
    "logP"
};
vector<string> DESCRIPTOR_API physicochemical_names(str_physicochemical, str_physicochemical+1);

// names for charged partial surface area (CPSA)
// different kind of atomic charge could yield various such descriptors
static string str_cpsa[] = {"TMSA (total molecular surface area)",
    "PPSA-1 (partial positive surface area)","PNSA-1 (partial positive surface area)",
    "DPSA-1 (PPSA-1 - PNSA-1)","FPSA-1 (PPSA-1 / TMSA)","FNSA-1 (PNSA-1 / TMSA)",
    "WPSA-1 (PPSA-1 * TMSA / 1000)","WNSA-1 (PNSA-1 * TMSA / 1000)",
    "PPSA-2 (total charge weighted PPSA)","PNSA-2 (total charge weighted PNSA)",
    "DPSA-2 (PPSA-2 - PNSA-2)","FPSA-2 (PPSA-2 / TMSA)","FNSA-2 (PNSA-2 / TMSA)",
    "WPSA-2 (PPSA-2 * TMSA / 1000)","WNSA-2 (PNSA-2 * TMSA / 1000)",
    "PPSA-3 (atomic charge weighted PPSA)","PNSA-3 (atomic charge weighted PNSA)",
    "DPSA-3 (PPSA-3 - PNSA-3)","FPSA-3 (PPSA-3 / TMSA)","FNSA-3 (PNSA-3 / TMSA)",
    "WPSA-3 (PPSA-3 * TMSA / 1000)","WNSA-3 (PNSA-3 * TMSA / 1000)",
    "RPCG (Qm.pos. / Qtot.pos.)","RPCS (SAm.pos * RPCG)",
    "RNCG (Qm.neg. / Qtot.neg.)","RNCS (SAm.neg * RPCG)",
    "HDSA (H-donors surface area)","FHDSA (HDSA / TMSA)",
    "HASA (H-acceptors surface area)","FHASA (HASA / TMSA)",
    "HBSA (HDSA + HASA)","FHBSA (HBSA / TMSA)",
    "HDCA (H-donors charge weighted surface area)","FHDCA (HDCA / TMSA)",
    "HACA (H-acceptors charge weighted surface area)","FHACA (HACA / TMSA)",
    "HBCA (HDCA + HACA)","FHBCA (HBCA / TMSA)"
};
vector<string> DESCRIPTOR_API cpsa_names(str_cpsa, str_cpsa+38);



string reference = "\
Todeschini R. and Consonni V. Molecular descriptors for chemoinformatics.\n\
  2nd Edition. Wiley\n\
CODESSA reference manual\n\
Mopac7, necessary if Quantum-Chemical descriptors are to be calculated\n\
  http://mopac7.sourceforge.net/\n\
OpenBabel, used to parse input molecules\n\
  http://sourceforge.net/projects/openbabel/\n\
  the version " + string(BABEL_VERSION) + " was used to compile this program.\n\
CPSA\n\
  Stanton D.T. and Jurs P.C. Anal. Chem. 1990, 62:2323-2329\n\
Balaban index\n\
  Ivanciuc O. et al. J. Chem. Inf. Comput. Sci. 1998, 38:395-401\n\
Randic index\n\
  Randic M. J. Am. Chem. Soc. 1975, 97:6609-6615\n\
  Randic M. J. Chem. Inf. Comput. Sci. 2001, 41:607-613\n\
Sanderson electronegativity\n\
  Sanderson R.T. J. Am. Chem. Soc. 1983. 105:2259-2261\n\
TPSA\n\
  Ertl P., Rohde B., and Selzer P. J. Med. Chem. 2000. 43:3714-3717\n\
Zefirov's PC\n\
  Oliferenko A.A. et al. J. Phys. Org. Chem. 2001. 14:355-369\n\
Quantum-Chemical descriptors\n\
  Karelson M. and Lobanov V.S. Chem. Rev. 1996. 96:1027-1043\n\
  Simas A.M. Theoret. Chim. Acta. 1982. 62:1-16\n\
Mulliken population analysis\n\
  Mulliken R.S. J. Chem. Phys. 1955. 10: 1833-1840, 1841-1846\n\
  Mulliken R.S. J. Chem. Phys. 1955. 12: 2338-2342, 2343-2346\n\
Shadow\n\
  Stouch T.R. and Jurs P.C. J. Chem. Inf. Comput. Sci. 1986. 26:4-12\n\
molecular volume\n\
  Richards F.M. Ann. Rev. Biophys. Bioeng. 1977. 6:151-176\n\
solvent accessible surface area\n\
  Shrake A. and Rupley J.A. J. Mol. Biol. 1973. 79:351-371\n\
E-state index\n\
  Hall L.H. and Kier L.B. J. Chem. Inf. Comput. Sci. 1995. 35:1039-1045\n\
logP\n\
  Wildman, S.A. and Crippen, G.M., J. Chem. Inf. Comput. Sci., 1999. 39:868-873\n\
Parameters\n\
  http://en.wikipedia.org/wiki/Electronegativity\n\
  R. T. Sanderson. J. Am. Chem. Soc. 1983, 105, 2259-2261\n\
  Cordero B. et al. Dalton. Trans., 2008: 2832-2838\n\
";

