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
#include "PAKViscosite.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

PAKViscosite::PAKViscosite(){}

//***********************************************************************

PAKViscosite::PAKViscosite(int& nombreGPA, Eos** eos, int &nombrePhases, string nomFichier)
{
  m_muk = new double[nombrePhases];
  for (int k = 0; k < nombrePhases; k++) {
    m_muk[k] = eos[k]->getMu();
  }
  m_numGPA = nombreGPA++;
}

//***********************************************************************

PAKViscosite::~PAKViscosite(){ delete[] m_muk; }


//***********************************************************************

void PAKViscosite::ajouteGrandeurPhysAdd(Cellule *cell)
{
  cell->getVecGrandeursPhysAdd().push_back(new GPAViscosite(this));
}

//***********************************************************************

void PAKViscosite::resolFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases)
{
  Cellule *cellGauche;
  Cellule *cellDroite;
  cellGauche = bord->getCellGauche();
  cellDroite = bord->getCellDroite();

  // Recopie vitesses et gradients des cellules gauche et droite
  m_vitesseGauche = cellGauche->getMelange()->getVitesse();
  m_vitesseDroite = cellDroite->getMelange()->getVitesse();

  m_gradUGauche = cellGauche->getGPA(m_numGPA)->getGrad(1);
  m_gradUDroite = cellDroite->getGPA(m_numGPA)->getGrad(1);
  m_gradVGauche = cellGauche->getGPA(m_numGPA)->getGrad(2);
  m_gradVDroite = cellDroite->getGPA(m_numGPA)->getGrad(2);
  m_gradWGauche = cellGauche->getGPA(m_numGPA)->getGrad(3);
  m_gradWDroite = cellDroite->getGPA(m_numGPA)->getGrad(3);

  // Calcul du mu de melange a gauche et a droite
  double muMelGauche(0.), muMelDroite(0.);
  for (int k = 0; k < nombrePhases; k++) {
    muMelGauche += cellGauche->getPhase(k)->getAlpha()*m_muk[k];
    muMelDroite += cellDroite->getPhase(k)->getAlpha()*m_muk[k];
  }

  Face *face;
  face = bord->getFace();
  m_normale = face->getNormale();
  m_tangente = face->getTangente();
  m_binormale = face->getBinormale();

  //Projection des vitesses et gradients sur repere attache a la face
  m_vitesseGauche.projection(m_normale, m_tangente, m_binormale);
  m_vitesseDroite.projection(m_normale, m_tangente, m_binormale);
  m_gradUGauche.projection(m_normale, m_tangente, m_binormale);
  m_gradUDroite.projection(m_normale, m_tangente, m_binormale);
  m_gradVGauche.projection(m_normale, m_tangente, m_binormale);
  m_gradVDroite.projection(m_normale, m_tangente, m_binormale);
  m_gradWGauche.projection(m_normale, m_tangente, m_binormale);
  m_gradWDroite.projection(m_normale, m_tangente, m_binormale);

  //Calcul des termes supplementaires du gradient dans la direction de la normale a la face.
  //BordDeMaille *bordCell(0);
  //Cellule *A1(0), *A3(0), *B1(0), *B3(0);
  //for (int b = 0; b < cellGauche->getBordsSize(); b++) {
  //  bordCell = cellGauche->getBord(b);
  //  string typeBord = bordCell->quiSuisJe();
  //  if (typeBord == "INTERNE")
  //  {
  //    if (fabs(bordCell->getFace()->getNormale().scalaire(m_normale)) < 1.e-6) {
  //      if (bordCell->getCellDroite() != cellGauche) { A1 = bordCell->getCellDroite(); }
  //      else if (bordCell->getCellGauche() != cellGauche) { A3 = bordCell->getCellGauche(); }
  //    }
  //  }
  //}
  //for (int b = 0; b < cellDroite->getBordsSize(); b++) {
  //  bordCell = cellDroite->getBord(b);
  //  string typeBord = bordCell->quiSuisJe();
  //  if (typeBord == "INTERNE")
  //  {
  //    if (fabs(bordCell->getFace()->getNormale().scalaire(m_normale)) < 1.e-6) {
  //      if (bordCell->getCellDroite() != cellDroite) { B1 = bordCell->getCellDroite(); }
  //      else if (bordCell->getCellGauche() != cellDroite) { B3 = bordCell->getCellGauche(); }
  //    }
  //  }
  //}

  //double termeSup = 0.;
  //if (A1 != 0 && B1 != 0) {
  //  termeSup = (B1->getMelange()->getVitesse().getX() - A1->getMelange()->getVitesse().getX()) / A1->distance(B1);
  //  m_gradUGauche.setX(m_gradUGauche.getX() + termeSup);
  //  m_gradUDroite.setX(m_gradUDroite.getX() + termeSup);

  //  termeSup = (B1->getMelange()->getVitesse().getY() - A1->getMelange()->getVitesse().getY()) / A1->distance(B1);
  //  m_gradVGauche.setX(m_gradVGauche.getX() + termeSup);
  //  m_gradVDroite.setX(m_gradVDroite.getX() + termeSup);
  //}
  //if (A3 != 0 && B3 != 0) {
  //  termeSup = (B3->getMelange()->getVitesse().getX() - A3->getMelange()->getVitesse().getX()) / A3->distance(B3);
  //  m_gradUGauche.setX(m_gradUGauche.getX() + termeSup);
  //  m_gradUDroite.setX(m_gradUDroite.getX() + termeSup);

  //  termeSup = (B3->getMelange()->getVitesse().getY() - A3->getMelange()->getVitesse().getY()) / A3->distance(B3);
  //  m_gradVGauche.setX(m_gradVGauche.getX() + termeSup);
  //  m_gradVDroite.setX(m_gradVDroite.getX() + termeSup);
  //}

  //Distances cellules/bord pour ponderation sur le flux
  double distGauche = cellGauche->distance(bord);
  double distDroite = cellDroite->distance(bord);

  this->resolFluxViscositeInterne(m_vitesseGauche, m_vitesseDroite, m_gradUGauche, m_gradUDroite, 
    m_gradVGauche, m_gradVDroite, m_gradWGauche, m_gradWDroite, muMelGauche, muMelDroite, distGauche, distDroite, nombrePhases);

  //Projection du flux sur le repere absolu
  bord->getMod()->projectionRepereAbsolu(m_normale, m_tangente, m_binormale);
}

//***********************************************************************

void PAKViscosite::resolFluxPhysAddLimite(BordDeMaille *bord, const int &nombrePhases)
{
  ////KS//DEV//VISC CL Inj, Res, Sortie a faire

  Cellule *cellGauche;
  cellGauche = bord->getCellGauche();

  //Recopie vitesses et gradients des cellules gauche et droite
  m_vitesseGauche = cellGauche->getMelange()->getVitesse();
  m_gradUGauche = cellGauche->getGPA(m_numGPA)->getGrad(1);
  m_gradVGauche = cellGauche->getGPA(m_numGPA)->getGrad(2);
  m_gradWGauche = cellGauche->getGPA(m_numGPA)->getGrad(3);

  // Calcul du mu de melange a gauche et a droite
  double muMelGauche(0.);
  for (int k = 0; k < nombrePhases; k++) {
    muMelGauche += cellGauche->getPhase(k)->getAlpha()*m_muk[k];
  }

  Face *face;
  face = bord->getFace();
  m_normale = face->getNormale();
  m_tangente = face->getTangente();
  m_binormale = face->getBinormale();

  //Projection des vitesses et gradients sur repere attache a la face
  m_vitesseGauche.projection(m_normale, m_tangente, m_binormale);
  m_gradUGauche.projection(m_normale, m_tangente, m_binormale);
  m_gradVGauche.projection(m_normale, m_tangente, m_binormale);
  m_gradWGauche.projection(m_normale, m_tangente, m_binormale);

  //Distances cellules/bord pour ponderation sur le flux
  double distGauche = cellGauche->distance(bord);

  int typeBord = bord->quiSuisJe();
  if (typeBord == 1) { this->resolFluxViscositeAbs(m_vitesseGauche, m_gradUGauche, m_gradVGauche, m_gradWGauche, muMelGauche, distGauche, nombrePhases); }
  else if (typeBord == 2) { this->resolFluxViscositeMur(m_vitesseGauche, muMelGauche, nombrePhases, distGauche); }
  else if (typeBord == 3) { this->resolFluxViscositeAbs(m_vitesseGauche, m_gradUGauche, m_gradVGauche, m_gradWGauche, muMelGauche, distGauche, nombrePhases); }
  else if (typeBord == 4) { this->resolFluxViscositeAbs(m_vitesseGauche, m_gradUGauche, m_gradVGauche, m_gradWGauche, muMelGauche, distGauche, nombrePhases); }
  else { this->resolFluxViscositeAutres(m_vitesseGauche, m_gradUGauche, m_gradVGauche, m_gradWGauche, muMelGauche, distGauche, nombrePhases); }
  // etc... CL pas gerees pour la capillarite, faire attention

  //Projection du flux sur le repere absolu
  bord->getMod()->projectionRepereAbsolu(m_normale, m_tangente, m_binormale);
}

//***********************************************************************

void PAKViscosite::resolFluxViscositeInterne(Coord &vitesseGauche, Coord &vitesseDroite, Coord &gradUGauche, Coord &gradUDroite, 
  Coord &gradVGauche, Coord &gradVDroite, Coord &gradWGauche, Coord &gradWDroite, double &muMelGauche, double &muMelDroite, double &distGauche, double &distDroite, int nombrePhases) const
{
	//Extraction des donnees
	double uL, vL, wL, uR, vR, wR;
	double du1L, du2L, du3L, du1R, du2R, du3R;
  double dv1L, dv2L, dv1R, dv2R;
  double dw1L, dw3L, dw1R, dw3R;
	uL = vitesseGauche.getX();
	vL = vitesseGauche.getY();
	wL = vitesseGauche.getZ();
	uR = vitesseDroite.getX();
	vR = vitesseDroite.getY();
	wR = vitesseDroite.getZ();

	du1L = gradUGauche.getX();
  du2L = gradUGauche.getY();
  du3L = gradUGauche.getZ();
  du1R = gradUDroite.getX();
  du2R = gradUDroite.getY();
  du3R = gradUDroite.getZ();

  dv1L = gradVGauche.getX();
  dv2L = gradVGauche.getY();
  dv1R = gradVDroite.getX();
  dv2R = gradVDroite.getY();

  dw1L = gradWGauche.getX();
  dw3L = gradWGauche.getZ();
  dw1R = gradWDroite.getX();
  dw3R = gradWDroite.getZ();

	//Donnees du bord de maille
  double u, v, w;
  double du1, du2, du3;
  double dv1, dv2;
  double dw1, dw3;
  double muMel;
	u = (uL*distDroite + uR*distGauche) / (distGauche + distDroite);
	v = (vL*distDroite + vR*distGauche) / (distGauche + distDroite);
	w = (wL*distDroite + wR*distGauche) / (distGauche + distDroite);

	du1 = (du1L*distDroite + du1R*distGauche) / (distGauche + distDroite);
  du2 = (du2L*distDroite + du2R*distGauche) / (distGauche + distDroite);
  du3 = (du3L*distDroite + du3R*distGauche) / (distGauche + distDroite);

  dv1 = (dv1L*distDroite + dv1R*distGauche) / (distGauche + distDroite);
  dv2 = (dv2L*distDroite + dv2R*distGauche) / (distGauche + distDroite);

  dw1 = (dw1L*distDroite + dw1R*distGauche) / (distGauche + distDroite);
  dw3 = (dw3L*distDroite + dw3R*distGauche) / (distGauche + distDroite);

  muMel = (muMelGauche*distDroite + muMelDroite*distGauche) / (distGauche + distDroite); //FP//Q// Pas compris ca...

	//Ecriture des termes visqueux sur chacunes des equations de fluxTempXXX
	for (int k = 0; k<nombrePhases; k++)
	{
	  fluxTempKapila->m_alpha[k] = 0.;
	  fluxTempKapila->m_masse[k] = 0.;
	  fluxTempKapila->m_energ[k] = 0.;
	}
  fluxTempKapila->m_qdm.setX(-muMel / 3. * (4.*du1 - 2.*(dv2 + dw3)));
  fluxTempKapila->m_qdm.setY(-muMel * (dv1 + du2));
  fluxTempKapila->m_qdm.setZ(-muMel * (dw1 + du3));
  fluxTempKapila->m_energMelange = -muMel * (4. / 3.*du1*u + (dv1 + du2)*v + (dw1 + du3)*w - 2. / 3.*(dv2 + dw3)*u);
  //fluxTempKapila->m_qdm.setX(0.);
  //fluxTempKapila->m_qdm.setY(0.);
  //fluxTempKapila->m_qdm.setZ(0.);
  //fluxTempKapila->m_energMelange = 0.;
}

//***********************************************************************

void PAKViscosite::resolFluxViscositeAbs(Coord &vitesseGauche, Coord &gradUGauche, Coord &gradVGauche, Coord &gradWGauche, double &muMelGauche, double &distGauche, int nombrePhases) const
{
  this->resolFluxViscositeInterne(vitesseGauche, vitesseGauche, gradUGauche, gradUGauche,
    gradVGauche, gradVGauche, gradWGauche, gradWGauche, muMelGauche, muMelGauche, distGauche, distGauche, nombrePhases);
}

//***********************************************************************

void PAKViscosite::resolFluxViscositeMur(Coord &vitesseGauche, double &muMelGauche, int nombrePhases, double &distGauche) const
{
  double du1 = -vitesseGauche.getX() / distGauche;
  double dv1 = -vitesseGauche.getY() / distGauche;
  double dw1 = -vitesseGauche.getZ() / distGauche;

  for (int k = 0; k<nombrePhases; k++)
  {
    fluxTempKapila->m_alpha[k] = 0.;
    fluxTempKapila->m_masse[k] = 0.;
    fluxTempKapila->m_energ[k] = 0.;
  }
  fluxTempKapila->m_qdm.setX(-muMelGauche / 3. * 4. * du1 );
  fluxTempKapila->m_qdm.setY(-muMelGauche * dv1);
  fluxTempKapila->m_qdm.setZ(-muMelGauche * dw1);
  fluxTempKapila->m_energMelange = 0.;
}

//***********************************************************************

void PAKViscosite::resolFluxViscositeAutres(Coord &vitesseGauche, Coord &gradUGauche, Coord &gradVGauche, Coord &gradWGauche, double &muMelGauche, double &distGauche, int nombrePhases) const
{
  //Pas gerer pour le moment
  cout << "CL visqueuse non geree" << endl;

  // Pour eviter bug quand non gerer
  for (int k = 0; k<nombrePhases; k++) {
    fluxTempKapila->m_alpha[k] = 0.;
    fluxTempKapila->m_masse[k] = 0.;
    fluxTempKapila->m_energ[k] = 0.;
  }
  fluxTempKapila->m_qdm = 0.;
  fluxTempKapila->m_energMelange = 0.;
}

//***********************************************************************

void PAKViscosite::ajoutNonCons(Cellule *cell, const int &nombrePhases)
{
  double du1 = cell->getGPA(m_numGPA)->getGrad(1).getX();
  double dv2 = cell->getGPA(m_numGPA)->getGrad(2).getY();
  double dw3 = cell->getGPA(m_numGPA)->getGrad(3).getZ();
  double termeNonCons = 2.*(du1*du1 + dv2*dv2 + dw3*dw3) - 2. / 3.*(du1 + dv2 + dw3)*(du1 + dv2 + dw3);

  for (int k = 0; k<nombrePhases; k++) {
    fluxTempKapila->m_alpha[k] = 0.;
    fluxTempKapila->m_masse[k] = 0.;
    fluxTempKapila->m_energ[k] = cell->getPhase(k)->getAlpha()*m_muk[k] * termeNonCons;
  }
  fluxTempKapila->m_qdm = 0.;
  fluxTempKapila->m_energMelange = 0.;

  cell->getCons()->ajoutFlux(1., nombrePhases);
}

//***********************************************************************

void PAKViscosite::communicationsPhysAdd(Cellule **cellules, const int &dim)
{
  Calcul_Parallele.communicationsVecteur(cellules, "GPA", dim, m_numGPA, 1); //m_gradU
  Calcul_Parallele.communicationsVecteur(cellules, "GPA", dim, m_numGPA, 2); //m_gradV
  Calcul_Parallele.communicationsVecteur(cellules, "GPA", dim, m_numGPA, 3); //m_gradW
}

//***********************************************************************

void PAKViscosite::communicationsPhysAddAMR(Cellule **cellules, const int &dim, const int &lvl)
{
	Calcul_Parallele.communicationsVecteurAMR(cellules, "GPA", dim, lvl, m_numGPA, 1); //m_gradU
	Calcul_Parallele.communicationsVecteurAMR(cellules, "GPA", dim, lvl, m_numGPA, 2); //m_gradV
	Calcul_Parallele.communicationsVecteurAMR(cellules, "GPA", dim, lvl, m_numGPA, 3); //m_gradW
}

//***********************************************************************
