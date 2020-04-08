#include "Objet_math.h"

void SELEM1(VECDBL &vdle, MATDBL &vcore, VECDBL &vmate, MATDBL &sige);

void SELEM2(VECDBL &vdle, MATDBL &vcore, VECDBL &vmate, MATDBL &sige);

//
//-----------------------------------------------------------------------
void gradient(int type, int nelt, int ndle, int ndln, int nnel, int ndim, int nsig, int nphy,
              VECDBL &ddl, MATDBL &vmat, VECINT &kmat, MATDBL &coord,
              MATINT &nelem, TRIDBL &sigma) {
//=======================================================================
//  calcul des contraintes aux noeuds selon le type d'éléments
//-----------------------------------------------------------------------
    int i, j, k, n, ieq, m, km, iel;
    VECDBL vdle(ndle), vmate(nphy);
    MATDBL sige(nsig, nnel), vcore(nnel, ndim);
    for (iel = 0; iel < nelt; iel++) {
        km = kmat[iel];
        for (i = 0; i < nphy; i++) vmate[i] = vmat[km - 1][i];
        m = 0;
        for (j = 0; j < nnel; j++) {
            n = nelem[iel][j];
            for (k = 0; k < ndim; k++) vcore[j][k] = coord[n - 1][k];
            for (k = 0; k < ndln; k++) {
                ieq = (n - 1) * ndln + k;
                vdle[m++] = ddl[ieq];
            }
        }
        switch (type) {
            case 1  :
                SELEM1(vdle, vcore, vmate, sige);
                break;
            case 2  :
                SELEM2(vdle, vcore, vmate, sige);
                break;
            default :
                cout << "L'élément du type  " << type << "  n'est pas implanté dans gradient" << endl;
        }
        for (k = 0; k < nsig; k++) for (j = 0; j < nnel; j++) sigma(k, j, iel) = sige[k][j];
    }
}

void SELEM1(VECDBL &vdle, MATDBL &vcore, VECDBL &vmate, MATDBL &sige) {
//=======================================================================
//     Calcul des contraintes aux noeuds
//     Elément xxxxx
//-----------------------------------------------------------------------

}

void SELEM2(VECDBL &vdle, MATDBL &vcore, VECDBL &vmate, MATDBL &sige) {
//=======================================================================
//     Calcul des contraintes aux noeuds
//     Elément xxxxx
//-----------------------------------------------------------------------

}
