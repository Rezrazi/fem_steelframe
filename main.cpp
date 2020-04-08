#include <iomanip>
#include <cstring>
#include "Objet_math.h"

using namespace std;

void elemlib(int iel, int type, MATDBL &vcore, VECDBL &vmate, MATDBL &vke, VECDBL &vfe);

void gradient(int type, int nelt, int ndle, int ndln, int nnel, int ndim, int nsig, int nphy,
              VECDBL &ddl, MATDBL &vmat, VECINT &kmat, MATDBL &coord,
              MATINT &nelem, TRIDBL &sigma);

void assemblage(VECINT &loce, MATDBL &vke, VECDBL &vfe, MATDBL &vkg, VECDBL &vfg);

void modifKF(MATINT &indice, MATDBL &valimp, MATDBL &vkg, VECDBL &vfg);

void identifier(int type, fstream &Out, int &idasc, char *dofname,
                int &ndim, int &ndln, int &nnel, int &nsig) {
    switch (type) {
        case 1  :
            Out << "\tType d'éléments : T3 " << endl;
            idasc = 23;
            ndim = 2;
            ndln = 2;
            nnel = 3;
            nsig = 4;
            strcpy_s(dofname, 100, "6  UX   UY   Sxx  Syy  Sxy Mises");
            break;
        case 2  :
            Out << "\tType d'éléments : Q4 " << endl;
            idasc = 24;
            ndim = 2;
            ndln = 2;
            nnel = 4;
            nsig = 4;
            strcpy_s(dofname, 100, "6  UX   UY   Sxx  Syy  Sxy Mises");
            break;
        default :
            Out << "L'element type  " << type << "  n'est pas encore implante" << endl;
            exit(0);
    }
}

int main() {
    char Nom_du_probleme[30], fich_inp[30], fich_out[30], fich_asc[30];
    cout << "nom du fichier svp? ";
    cin >> Nom_du_probleme;
    strcpy_s(fich_inp, Nom_du_probleme);
    strcat_s(fich_inp, ".inp");
    strcpy_s(fich_out, Nom_du_probleme);
    strcat_s(fich_out, ".out");
    strcpy_s(fich_asc, Nom_du_probleme);
    strcat_s(fich_asc, ".asc");

    //------ ouvrire le fichier de donnees et creer deux fichiers de resultats

    fstream Inp, Out, Asc;
    Inp.open(fich_inp, ios::in);
    if (!Inp) {
        cout << "erreur d ouverture du fichier .inp\n";
        exit(0);
    }
    Out.open(fich_out, ios::out);
    if (!Out) {
        cout << "erreur d ouverture du fichier .out\n";
        exit(0);
    }
    Asc.open(fich_asc, ios::out);
    if (!Asc) {
        cout << "erreur d ouverture du fichier .asc\n";
        exit(0);
    }

    Out << "   ================================================================\n";
    Out << "  |                                                                |\n";
    Out << "  |                             CCMEF                              |\n";
    Out << "  |         Code de Calcul par la Méthode des Eléments Finis       |\n";
    Out << "  |                                                                |\n";
    Out << "  |            Réalise par :                                       |\n";
    Out << "  |                   Date :                                       |\n";
    Out << "  |                                                                |\n";
    Out << "  |                                                                |\n";
    Out << "  |                   Enseignant : Zhi-Qiang FENG                  |\n";
    Out << "  |                                                                |\n";
    Out << "   ================================================================\n\n\n";

    int i, j, k, h, id, idasc, n, km, iel, nnt, nelt, nnel, ndln, ndim, nmat, nphy, nsig, ndle, ieq, neq, type;
    double force;
    char delim = '\n', ligne[200], titre[200], dofname[100];
    Inp.getline(titre, 200, delim);
    Out << "============ " << titre << endl << endl;

    Inp.getline(ligne, 200, delim);
    Inp >> type;
    identifier(type, Out, idasc, dofname, ndim, ndln, nnel, nsig);
    Out << endl;

    Inp >> ligne;
    Out << "============ " << ligne << endl;

    //-..... paramètres du problème

    Inp >> nnt >> nelt >> nmat >> nphy;
    ndle = nnel * ndln;
    neq = nnt * ndln;
    Out << endl;
    Out << "\tnombre total de noeuds:               nnt  = " << nnt << endl;
    Out << "\tdimension du problème:                ndim = " << ndim << endl;
    Out << "\tnombre total d'éléments:              nelt = " << nelt << endl;
    Out << "\tnombre de noeuds par élément:         nnel = " << nnel << endl;
    Out << "\tnombre ddls par noeud:                ndln = " << ndln << endl;
    Out << "\tnombre ddls par élément:              ndle = " << ndle << endl;
    Out << "\tnombre de matériaux:                  nmat = " << nmat << endl;
    Out << "\tnombre de propriétés physiques:       nphy = " << nphy << endl;
    Out << "\tnombre de composantes gradients:      nsig = " << nsig << endl;
    Out << "\tnombre d equations a résoudre:         neq = " << neq << endl;
    Out << endl;

    //------  allocation des tables
    MATDBL valimp(nnt, ndln), vmat(nmat, nphy), coord(nnt, ndim), vcore(nnel, ndim);
    TRIDBL sigma(nsig, nnel, nelt);
    MATINT indice(nnt, ndln), nelem(nelt, nnel);
    VECINT kmat(nelt), loce(ndle);
    VECDBL vfe(ndle), vmate(nphy);

    MATDBL vke(ndle, ndle), vkg(neq, neq);
    VECDBL vfg(neq), ddl(neq);

    //------	coordonnées des noeuds
    Inp >> ligne;
    Out << "============ " << ligne << " ============ " << endl;
    for (i = 0; i < nnt; i++) {
        Inp >> n;
        for (j = 0; j < ndim; j++) Inp >> coord[n - 1][j];
        Out << setw(5) << n;
        for (j = 0; j < ndim; j++) Out << setw(16) << coord[n - 1][j];
        Out << endl;
    }

    //------	connectivités des éléments
    Inp >> ligne;
    Out << "============ " << ligne << " ============ " << endl;
    for (iel = 0; iel < nelt; iel++) {
        Inp >> n;
        for (j = 0; j < nnel; j++) Inp >> nelem[n - 1][j];
        Inp >> kmat[n - 1];
        Out << setw(5) << n;
        for (j = 0; j < nnel; j++) Out << setw(6) << nelem[n - 1][j];
        Out << setw(6) << kmat[n - 1];
        Out << endl;
    }

    //------	propriétés matérielles et géométriques
    Inp >> ligne;
    Out << "============ " << ligne << " ============ " << endl;
    for (i = 0; i < nmat; i++) {
        Inp >> n;
        for (j = 0; j < nphy; j++) Inp >> vmat[n - 1][j];
        Out << setw(5) << n;
        for (j = 0; j < nphy; j++) Out << setw(16) << vmat[n - 1][j];
        Out << endl;
    }

    //------	condition aux limites
    Inp >> ligne;
    Out << "============ " << ligne << " ============ " << endl;
    int ndlt = 0;
    for (i = 0; i <= nnt; i++) {
        Inp >> n;
        if (n == -9999) break;
        for (j = 0; j < ndln; j++) {
            Inp >> indice[n - 1][j];
            if (indice[n - 1][j] == 1) ndlt++;
        }
        for (j = 0; j < ndln; j++) Inp >> valimp[n - 1][j];
        Out << setw(5) << n;
        for (j = 0; j < ndln; j++) Out << setw(6) << indice[n - 1][j];
        for (j = 0; j < ndln; j++) Out << setw(16) << valimp[n - 1][j];
        Out << endl;
    }
    Out << "\tnombre total des ddl imposés:         ndlt = " << ndlt << endl;
    Out << endl;

    //------    sollicitation concentrée aux noeuds
    Inp >> ligne;
    Out << "============ " << ligne << " ============ " << endl;
    for (i = 0; i < nnt; i++) {
        Inp >> n;
        if (n == -9999) break;
        Out << setw(5) << n;
        for (j = 0; j < ndln; j++) {
            ieq = (n - 1) * ndln + j;
            Inp >> vfg[ieq];
            Out << setw(16) << vfg[ieq];
        }
        Out << endl;
    }

    //==========Construction de K et F =====================
    //------	boucle sur les éléments
    for (iel = 0; iel < nelt; iel++) {
        km = kmat[iel];
        for (i = 0; i < nphy; i++) vmate[i] = vmat[km - 1][i];
        h = 0;
        for (j = 0; j < nnel; j++) {
            n = nelem[iel][j];
            id = (n - 1) * ndln;
            for (k = 0; k < ndln; k++) loce[h++] = id + k;
            for (k = 0; k < ndim; k++) vcore[j][k] = coord[n - 1][k];
        }
        //------ calcul des matrices et vecteurs élémentaires
        elemlib(iel, type, vcore, vmate, vke, vfe);
        //------ assemblage
        assemblage(loce, vke, vfe, vkg, vfg);
    }

    //------ modification de K et F pour tenir compte des ddl.imposees

    modifKF(indice, valimp, vkg, vfg);

    //------ résolution du système [K]{u} = {F}
    vkg.solver_gauss(vfg, ddl);
    //------ impression des résultats
    Out << "\n\n======= Degrés de libert?aux noeuds =======\n";
    k = 0;
    for (i = 0; i < nnt; i++) {
        Out << endl << setw(5) << i + 1;
        for (j = 0; j < ndln; j++) {
            ieq = i * ndln + j;
            Out << setw(16) << ddl[ieq];
        }
    }

    //------	calcul des gradients
    Out << "\n\n======= Contraintes dans les éléments =======\n";
    gradient(type, nelt, ndle, ndln, nnel, ndim, nsig, nphy, ddl, vmat, kmat, coord, nelem, sigma);
    for (iel = 0; iel < nelt; iel++) {
        Out << " ------------------------ Elément : " << setw(5) << iel + 1 << endl;
        for (j = 0; j < nnel; j++) {
            n = nelem[iel][j];
            Out << setw(5) << n;
            for (k = 0; k < nsig; k++) Out << setw(16) << sigma(k, j, iel);
            Out << endl;
        }
    }
    //------    calcul des valeurs moyennes aux noeuds
    int ino;
    VECDBL valeur(nnt);
    VECINT nodedup(nnt);
    MATDBL sigma_node(nnt, nsig);
    for (k = 0; k < nsig; k++) {
        for (ino = 0; ino < nnt; ino++) {
            valeur[ino] = 0.;
            nodedup[ino] = 0;
        }
        for (iel = 0; iel < nelt; iel++) {
            for (j = 0; j < nnel; j++) {
                n = nelem[iel][j];
                valeur[n - 1] += sigma(k, j, iel);
                nodedup[n - 1]++;
            }
        }
        for (ino = 0; ino < nnt; ino++) {
            if (nodedup[ino])valeur[ino] /= nodedup[ino]; else valeur[ino] = 0.;
            sigma_node[ino][k] = valeur[ino];
        }
    }
    //====================================================================
    //  maillage dans le fichier .asc pour le postprocesseur FER/View
    //====================================================================
    Asc << titre << endl;
    Asc << nnt << "\t" << ndim << endl;
    for (i = 0; i < nnt; i++) {
        Asc << setw(5) << i + 1;
        for (j = 0; j < ndim; j++) Asc << setw(16) << coord[i][j] << "\t";
        Asc << endl;
    }
    Asc << nelt << endl;
    int nun = 1;
    for (iel = 0; iel < nelt; iel++) {
        km = kmat[iel];
        Asc << setw(5) << iel + 1 << setw(5) << nnel << setw(5) << idasc << setw(5) << km << setw(5) << nun;
        for (j = 0; j < nnel; j++) Asc << setw(5) << nelem[iel][j];
        Asc << endl;
    }
    //====================================================================
    //  résultats dans le fichier .asc pour le postprocesseur FER/View
    //====================================================================
    Asc << dofname << endl;
    Asc << "  1" << endl;
    Asc << "  1   0.0000e+000" << endl;
    for (i = 0; i < nnt; i++) {
        Asc << setw(5) << i + 1;
        for (j = 0; j < ndln; j++) {
            ieq = i * ndln + j;
            Asc << setw(16) << ddl[ieq];
        }
        for (j = 0; j < nsig; j++) Asc << setw(16) << sigma_node[i][j];
        Asc << endl;
    }
    return 0;
}