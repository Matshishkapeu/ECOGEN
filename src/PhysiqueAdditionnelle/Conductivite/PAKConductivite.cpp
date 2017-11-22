//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <cmath>
#include <algorithm>
#include "PAKConductivite.h"

using namespace std;

//***********************************************************************

PAKConductivite::PAKConductivite() {}

//***********************************************************************

PAKConductivite::PAKConductivite(int& nombreGPA, Eos** eos, int &nombrePhases, string nomFichier)
{
  m_lambdak = new double[nombrePhases];
  for (int k = 0; k < nombrePhases; k++) {
    m_lambdak[k] = eos[k]->getLambda();
  }
  m_numGPA = nombreGPA++;
}

//***********************************************************************

PAKConductivite::~PAKConductivite(){ delete[] m_lambdak; }

//***********************************************************************

void PAKConductivite::ajouteGrandeurPhysAdd(Cellule *cell)
{
  cell->getVecGrandeursPhysAdd().push_back(new GPAConductivite(this,cell->getNombrePhases()));
}

//***********************************************************************

void PAKConductivite::resolFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases)
{
  Cellule *cellGauche;
  Cellule *cellDroite;
  cellGauche = bord->getCellGauche();
  cellDroite = bord->getCellDroite();

  Face *face;
  face = bord->getFace();
  m_normale = face->getNormale();
  m_tangente = face->getTangente();
  m_binormale = face->getBinormale();

  // Mise a zero du fluxTempKapila
  for (int k = 0; k<nombrePhases; k++) {
    fluxTempKapila->m_alpha[k] = 0.;
    fluxTempKapila->m_masse[k] = 0.;
    fluxTempKapila->m_energ[k] = 0.;
  }
  fluxTempKapila->m_qdm = 0.;
  fluxTempKapila->m_energMelange = 0.;

  for (int numPhase = 0; numPhase < nombrePhases; numPhase++) {
    // Recopie et projection sur repere attache a la face des gradients des cellules gauche et droite
    m_gradTkGauche = cellGauche->getGPA(m_numGPA)->getGrad(numPhase);
    m_gradTkDroite = cellDroite->getGPA(m_numGPA)->getGrad(numPhase);
    m_gradTkGauche.projection(m_normale, m_tangente, m_binormale);
    m_gradTkDroite.projection(m_normale, m_tangente, m_binormale);

    // Extraction des alphak
    double alphakGauche = cellGauche->getPhase(numPhase)->getAlpha();
    double alphakDroite = cellDroite->getPhase(numPhase)->getAlpha();

    this->resolFluxConductiviteInterne(m_gradTkGauche, m_gradTkDroite, alphakGauche, alphakDroite, numPhase);
  }

  //Projection du flux sur le repere absolu
  bord->getMod()->projectionRepereAbsolu(m_normale, m_tangente, m_binormale);
}

//***********************************************************************

void PAKConductivite::resolFluxPhysAddLimite(BordDeMaille *bord, const int &nombrePhases)
{
  //KS//DEV//COND On ne fait rien aux limites avec la conductivite pour le moment, a gerer un jour

  Cellule *cellGauche;
  cellGauche = bord->getCellGauche();

  Face *face;
  face = bord->getFace();
  m_normale = face->getNormale();
  m_tangente = face->getTangente();
  m_binormale = face->getBinormale();

  // Mise a zero du fluxTempKapila (permet de faire ensuite la somme des effets conductifs pour les differentes phases)
  for (int k = 0; k<nombrePhases; k++) {
    fluxTempKapila->m_alpha[k] = 0.;
    fluxTempKapila->m_masse[k] = 0.;
    fluxTempKapila->m_energ[k] = 0.;
  }
  fluxTempKapila->m_qdm = 0.;
  fluxTempKapila->m_energMelange = 0.;

  for (int numPhase = 0; numPhase < nombrePhases; numPhase++) {
    // Recopie et projection sur repere attache a la face des gradients des cellules gauche et droite
    m_gradTkGauche = cellGauche->getGPA(m_numGPA)->getGrad(numPhase);
    m_gradTkGauche.projection(m_normale, m_tangente, m_binormale);

    // Extraction des alphak
    double alphakGauche = cellGauche->getPhase(numPhase)->getAlpha();

    int typeBord = bord->quiSuisJe();
    if (typeBord == 1) { this->resolFluxConductiviteAbs(m_gradTkGauche, alphakGauche, numPhase); }
    else if (typeBord == 2) { this->resolFluxConductiviteMur(m_gradTkGauche, alphakGauche, numPhase); }
    else if (typeBord == 3) { this->resolFluxConductiviteSortie(m_gradTkGauche, alphakGauche, numPhase); }
    else if (typeBord == 4) { this->resolFluxConductiviteInjection(m_gradTkGauche, alphakGauche, numPhase); }
    else { this->resolFluxConductiviteAutres(m_gradTkGauche, alphakGauche, numPhase); }
    // etc... CL pas gerees pour la conductivite, faire attention
  }

  //Projection du flux sur le repere absolu
  bord->getMod()->projectionRepereAbsolu(m_normale, m_tangente, m_binormale);
}

//***********************************************************************

void PAKConductivite::resolFluxConductiviteInterne(Coord &gradTkGauche, Coord &gradTkDroite, double &alphakL, double &alphakR, int &numPhase) const
{
  //Donnees du bord de maille
  double dTkdx, alphak;
  dTkdx = (gradTkGauche.getX() + gradTkDroite.getX()) / 2.;
  alphak = (alphakL + alphakR) / 2.;

  //Ecriture des termes conductifs sur chacunes des equations de fluxTempXXX
  fluxTempKapila->m_energ[numPhase] += -alphak*m_lambdak[numPhase]* dTkdx;
  fluxTempKapila->m_energMelange += -alphak*m_lambdak[numPhase] * dTkdx;
}

//***********************************************************************

void PAKConductivite::resolFluxConductiviteAbs(Coord &gradTkGauche, double &alphakL, int &numPhase) const
{
  this->resolFluxConductiviteInterne(gradTkGauche, gradTkGauche, alphakL, alphakL, numPhase);
}

//***********************************************************************

void PAKConductivite::resolFluxConductiviteMur(Coord &gradTkGauche, double &alphakL, int &numPhase) const
{
  //Pas gerer pour le moment, juste un exemple
  //KS//DEV//COND// A faire pour couche limite thermique !!! ...

  // Pour éviter bug quand non gerer
  fluxTempKapila->m_qdm.setX(fluxTempKapila->m_qdm.getX() + 0.);
  fluxTempKapila->m_qdm.setY(fluxTempKapila->m_qdm.getY() + 0.);
  fluxTempKapila->m_qdm.setZ(fluxTempKapila->m_qdm.getZ() + 0.);
  fluxTempKapila->m_energMelange += 0.;
}

//***********************************************************************

void PAKConductivite::resolFluxConductiviteSortie(Coord &gradTkGauche, double &alphakL, int &numPhase) const
{
  //Pas gerer pour le moment, juste un exemple

  // Pour éviter bug quand non gerer
  fluxTempKapila->m_qdm.setX(fluxTempKapila->m_qdm.getX() + 0.);
  fluxTempKapila->m_qdm.setY(fluxTempKapila->m_qdm.getY() + 0.);
  fluxTempKapila->m_qdm.setZ(fluxTempKapila->m_qdm.getZ() + 0.);
  fluxTempKapila->m_energMelange += 0.;
}

//***********************************************************************

void PAKConductivite::resolFluxConductiviteInjection(Coord &gradTkGauche, double &alphakL, int &numPhase) const
{
  //Pas gerer pour le moment, juste un exemple

  // Pour éviter bug quand non gerer
  fluxTempKapila->m_qdm.setX(fluxTempKapila->m_qdm.getX() + 0.);
  fluxTempKapila->m_qdm.setY(fluxTempKapila->m_qdm.getY() + 0.);
  fluxTempKapila->m_qdm.setZ(fluxTempKapila->m_qdm.getZ() + 0.);
  fluxTempKapila->m_energMelange += 0.;
}

//***********************************************************************

void PAKConductivite::resolFluxConductiviteAutres(Coord &gradTkGauche, double &alphakL, int &numPhase) const
{
  //Pas gerer pour le moment
  cout << "CL conductivite non geree" << endl;

  // Pour éviter bug quand non gerer
  fluxTempKapila->m_qdm.setX(fluxTempKapila->m_qdm.getX() + 0.);
  fluxTempKapila->m_qdm.setY(fluxTempKapila->m_qdm.getY() + 0.);
  fluxTempKapila->m_qdm.setZ(fluxTempKapila->m_qdm.getZ() + 0.);
  fluxTempKapila->m_energMelange += 0.;
}

//***********************************************************************

void PAKConductivite::communicationsPhysAdd(Cellule **cellules, const int &dim)
{
  for (int k = 0; k < cellules[0]->getNombrePhases(); k++) {
    Calcul_Parallele.communicationsVecteur(cellules, "GPA", dim, m_numGPA, k);
  }
}

//***********************************************************************

void PAKConductivite::communicationsPhysAddAMR(Cellule **cellules, const int &dim, const int &lvl)
{
	for (int k = 0; k < cellules[0]->getNombrePhases(); k++) {
		Calcul_Parallele.communicationsVecteurAMR(cellules, "GPA", dim, lvl, m_numGPA, k);
	}
}

//***********************************************************************

double PAKConductivite::getLambdak(int &numPhase) const
{
  return m_lambdak[numPhase];
}

//***********************************************************************