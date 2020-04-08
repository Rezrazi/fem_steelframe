#include "Objet_math.h"

void assemblage(VECINT &loce, MATDBL &vke, VECDBL &vfe, MATDBL &vkg, VECDBL &vfg) {
/*
c=====================================================================
c
c     assemblage d'une matrice et/ou d'un vecteur elementaire
c       entrees
c          loce    vecteur de localisation de l'element
c          vke     matrice elementaire ke(pleine)
c          vfe     vecteur elementaire fe
c       sorties
c          vkg     matrice globale (matrice pleine)
c          vfg     vecteur sollicitations global
c=====================================================================
*/
    int i, j, ig, jg;
    int ndle = loce.n;
    for (i = 0; i < ndle; i++) {
        ig = loce[i];
        for (j = 0; j < ndle; j++) {
            jg = loce[j];
            vkg[ig][jg] += vke[i][j];
        }
        vfg[ig] += vfe[i];
    }
}

void modifKF(MATINT &indice, MATDBL &valimp, MATDBL &vkg, VECDBL &vfg) {
//======================================================================
// imposition des conditions aux limites
// remarques: utilisation de la methode du terme diagonal dominant
//======================================================================
    int i, j, ieq;
    int nnt = indice.n, ndln = indice.m;
    double alfa = 1.0e40;
    for (i = 0; i < nnt; i++) {
        for (j = 0; j < ndln; j++) {
            if (indice[i][j] == 1) {
                ieq = i * ndln + j;
                vkg[ieq][ieq] += alfa;
                vfg[ieq] = valimp[i][j] * alfa;
            }
        }
    }
}
