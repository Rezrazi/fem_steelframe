#include <iomanip>
#include "Objet_math.h"

using namespace std;

void ELEM1(int iel, MATDBL &vcore, VECDBL &vmate, MATDBL &vke, VECDBL &vfe);

void ELEM2(int iel, MATDBL &vcore, VECDBL &vmate, MATDBL &vke, VECDBL &vfe);

//
void elemlib(int iel, int type, MATDBL &vcore, VECDBL &vmate, MATDBL &vke, VECDBL &vfe) {
    switch (type) {
        case 1 :
            ELEM1(iel, vcore, vmate, vke, vfe);
            break;
        case 2 :
            ELEM2(iel, vcore, vmate, vke, vfe);
            break;
        default :
            cout << "L'element type  " << type << "  n'est pas encore implante" << endl;
            exit(0);
    }
}

void ELEM1(int iel, MATDBL &vcore, VECDBL &vmate, MATDBL &vke, VECDBL &vfe) {
//=======================================================================
//     calcul de la matrice de rigidité élémentaire
//     Elément xxxxxx
//-----------------------------------------------------------------------
}

void ELEM2(int iel, MATDBL &vcore, VECDBL &vmate, MATDBL &vke, VECDBL &vfe) {
//=======================================================================
//     calcul de la matrice de rigidité élémentaire
//     Elément xxxxxx
//-----------------------------------------------------------------------
}
