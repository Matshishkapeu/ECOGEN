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
#include "PAKCapillarite.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

PAKCapillarite::PAKCapillarite(){}

//***********************************************************************
/*!
*  Constructeur physiqueAdd a partir d une lecture au format XML
*  ex : <donneesCapillarite transport="couleur1" sigma="800."/>
*/
PAKCapillarite::PAKCapillarite(XMLElement *element, int& nombreGPA, vector<string> nomTransports, string nomFichier)
{
  XMLElement *sousElement(element->FirstChildElement("donneesCapillarite"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesCapillarite", nomFichier, __FILE__, __LINE__);
  //Recuperation des attributs
  //--------------------------
  XMLError erreur;
  //Nom de l equation de transport associee
  m_nomTransportAssocie = sousElement->Attribute("transport");
  if (m_nomTransportAssocie == "") throw ErreurXMLAttribut("transport", nomFichier, __FILE__, __LINE__);
  //On associe tout de suite avec le numero correspondant dans m_vecTransport des cellules
  unsigned int t(0);
  for (t = 0; t < nomTransports.size(); t++) {
    if (m_nomTransportAssocie == nomTransports[t]) { break; }
  }
  if (t != nomTransports.size()) { m_numTransportAssocie = t; }
  else { Erreurs::messageErreur("Equation de transport non trouve pour Capillarite"); }

  //Valeur tension de surface
  erreur = sousElement->QueryDoubleAttribute("sigma", &m_sigma);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("sigma", nomFichier, __FILE__, __LINE__);
  //Mise en correspondance avec les Grandeurs physiques additionelles
  m_numGPAGradC = nombreGPA++;
}

//***********************************************************************

PAKCapillarite::~PAKCapillarite(){}

//***********************************************************************

void PAKCapillarite::ajouteGrandeurPhysAdd(Cellule *cell)
{
  cell->getVecGrandeursPhysAdd().push_back(new GPACapillarite(this));
}

//***********************************************************************

double PAKCapillarite::calculEnergiePhysAdd(GrandeursPhysAdd* GPA)
{
  double energieCap = m_sigma*GPA->getGrad().norme(); //preparation energie capillaire volumique (sigma*gradC)
  return energieCap;
}

//***********************************************************************

void PAKCapillarite::resolFluxPhysAdd(BordDeMaille *bord, const int& nombrePhases)
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

  // Recopie et projection sur repere attache a la face des vitesses des cellules gauche et droite
  m_vitesseGauche = cellGauche->getMelange()->getVitesse();
  m_vitesseDroite = cellDroite->getMelange()->getVitesse();
  m_vitesseGauche.projection(m_normale, m_tangente, m_binormale);
  m_vitesseDroite.projection(m_normale, m_tangente, m_binormale);

  // Mise a zero du fluxTempKapila
  fluxTempKapila->miseAZero(nombrePhases);

  // Recopie et projection sur repere attache a la face des gradients des cellules gauche et droite
  m_gradCGauche = cellGauche->getGPA(m_numGPAGradC)->getGrad();
  m_gradCDroite = cellDroite->getGPA(m_numGPAGradC)->getGrad();
  m_gradCGauche.projection(m_normale, m_tangente, m_binormale);
  m_gradCDroite.projection(m_normale, m_tangente, m_binormale);

  this->resolFluxCapillariteInterne(m_vitesseGauche, m_vitesseDroite, m_gradCGauche, m_gradCDroite);

  //Projection du flux sur le repere absolu
  bord->getMod()->projectionRepereAbsolu(m_normale, m_tangente, m_binormale);
}

//***********************************************************************

void PAKCapillarite::resolFluxPhysAddLimite(BordDeMaille *bord, const int &nombrePhases)
{
  //KS//DEV//CAP On ne fait rien aux limites avec la capillarite pour le moment, a gerer un jour

  Cellule *cellGauche;
  cellGauche = bord->getCellGauche();

  Face *face;
  face = bord->getFace();
  m_normale = face->getNormale();
  m_tangente = face->getTangente();
  m_binormale = face->getBinormale();

  // Recopie et projection sur repere attache a la face des vitesses des cellules gauche et droite
  m_vitesseGauche = cellGauche->getMelange()->getVitesse();
  m_vitesseGauche.projection(m_normale, m_tangente, m_binormale);

  // Mise a zero du fluxTempKapila (permet de faire ensuite la somme des effets capillaires pour les differentes combinaisons de phases)
  fluxTempKapila->miseAZero(nombrePhases);
  
  // Recopie et projection sur repere attache a la face des gradients des cellules gauche et droite
  m_gradCGauche = cellGauche->getGPA(m_numGPAGradC)->getGrad();
  m_gradCGauche.projection(m_normale, m_tangente, m_binormale);

  int typeBord = bord->quiSuisJe();
  if (typeBord == 1) { this->resolFluxCapillariteAbs(m_vitesseGauche, m_gradCGauche); }
  else if (typeBord == 2) { this->resolFluxCapillariteMur(m_gradCGauche); }
  else if (typeBord == 3) { this->resolFluxCapillariteSortie(m_vitesseGauche, m_gradCGauche); }
  else if (typeBord == 4) { this->resolFluxCapillariteInjection(m_vitesseGauche, m_gradCGauche); }
  else { this->resolFluxCapillariteAutres(m_vitesseGauche, m_gradCGauche); }
  // etc... CL pas gerees pour la capillarite, faire attention

  //Projection du flux sur le repere absolu
  bord->getMod()->projectionRepereAbsolu(m_normale, m_tangente, m_binormale);
}

//***********************************************************************

void PAKCapillarite::resolFluxCapillariteInterne(Coord &vitesseGauche, Coord &vitesseDroite, Coord &gradCGauche, Coord &gradCDroite) const
{
	//Extraction des donnees
	double uL, vL, wL, uR, vR, wR;
	double w1L, w2L, w3L, normeWL, w1R, w2R, w3R, normeWR;
	uL = vitesseGauche.getX();
	vL = vitesseGauche.getY();
	wL = vitesseGauche.getZ();
	uR = vitesseDroite.getX();
	vR = vitesseDroite.getY();
	wR = vitesseDroite.getZ();
	w1L = gradCGauche.getX();
	w2L = gradCGauche.getY();
	w3L = gradCGauche.getZ();
	normeWL = gradCGauche.norme();
	w1R = gradCDroite.getX();
	w2R = gradCDroite.getY();
	w3R = gradCDroite.getZ();
	normeWR = gradCDroite.norme();

	//Donnees du bord de maille
	double u, v, w, w1, w2, w3, normeW;
	u = (uL + uR) / 2.;
	v = (vL + vR) / 2.;
	w = (wL + wR) / 2.;
	w1 = (w1L + w1R) / 2.;
	w2 = (w2L + w2R) / 2.;
	w3 = (w3L + w3R) / 2.;
	normeW = (normeWL + normeWR) / 2.;

	//Ecriture des termes capillaires sur chacunes des equations de fluxTempXXX
	if (normeW > 1.e-6) {
	  fluxTempKapila->m_qdm.setX(fluxTempKapila->m_qdm.getX() + m_sigma*(w1*w1 / normeW - normeW));
	  fluxTempKapila->m_qdm.setY(fluxTempKapila->m_qdm.getY() + m_sigma *w1*w2 / normeW);
	  fluxTempKapila->m_qdm.setZ(fluxTempKapila->m_qdm.getZ() + m_sigma *w1*w3 / normeW);
    fluxTempKapila->m_energMelange += m_sigma / normeW*(w1*w1*u + w1*w2*v + w1*w3*w);
	}
}

//***********************************************************************

void PAKCapillarite::resolFluxCapillariteAbs(Coord &vitesseGauche, Coord &gradCGauche) const
{
  this->resolFluxCapillariteInterne(vitesseGauche, vitesseGauche, gradCGauche, gradCGauche);
}

//***********************************************************************

void PAKCapillarite::resolFluxCapillariteMur(Coord &gradCGauche) const
{
  // Considere comme une symetrie et non comme un mur avec point triple capillaire

  //Donnees du bord de maille
  double normeW;
  normeW = gradCGauche.norme();

  //Ecriture des termes capillaires sur chacunes des equations de fluxTempXXX
  if (normeW > 1.e-6) {
    fluxTempKapila->m_qdm.setX(fluxTempKapila->m_qdm.getX() - m_sigma*normeW);
  }
}

//***********************************************************************

void PAKCapillarite::resolFluxCapillariteSortie(Coord &vitesseGauche, Coord &gradCGauche) const
{
  //Pas gerer pour le moment, juste un exemple

  // Pour éviter bug quand non gerer
  fluxTempKapila->m_qdm.setX(fluxTempKapila->m_qdm.getX() + 0.);
  fluxTempKapila->m_qdm.setY(fluxTempKapila->m_qdm.getY() + 0.);
  fluxTempKapila->m_qdm.setZ(fluxTempKapila->m_qdm.getZ() + 0.);
  fluxTempKapila->m_energMelange += 0.;
}

//***********************************************************************

void PAKCapillarite::resolFluxCapillariteInjection(Coord &vitesseGauche, Coord &gradCGauche) const
{
  //Pas gerer pour le moment, juste un exemple

  // Pour éviter bug quand non gerer
  fluxTempKapila->m_qdm.setX(fluxTempKapila->m_qdm.getX() + 0.);
  fluxTempKapila->m_qdm.setY(fluxTempKapila->m_qdm.getY() + 0.);
  fluxTempKapila->m_qdm.setZ(fluxTempKapila->m_qdm.getZ() + 0.);
  fluxTempKapila->m_energMelange += 0.;
}

//***********************************************************************

void PAKCapillarite::resolFluxCapillariteAutres(Coord &vitesseGauche, Coord &gradCGauche) const
{
  //Pas gerer pour le moment
  cout << "CL capillaire non geree" << endl;

  // Pour éviter bug quand non gerer
  fluxTempKapila->m_qdm.setX(fluxTempKapila->m_qdm.getX() + 0.);
  fluxTempKapila->m_qdm.setY(fluxTempKapila->m_qdm.getY() + 0.);
  fluxTempKapila->m_qdm.setZ(fluxTempKapila->m_qdm.getZ() + 0.);
  fluxTempKapila->m_energMelange += 0.;
}

//***********************************************************************

void PAKCapillarite::reinitialiseFonctionCouleur(vector<Cellule *> *cellulesLvl, int &lvl)
{
	for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) {
		if (!cellulesLvl[lvl][i]->getSplit()) { cellulesLvl[lvl][i]->reinitialiseFonctionCouleur(); }
	}
}

//***********************************************************************

void PAKCapillarite::communicationsPhysAdd(Cellule **cellules, const int &dim)
{
  Calcul_Parallele.communicationsVecteur(cellules, "GPA", dim, m_numGPAGradC);
}

//***********************************************************************

void PAKCapillarite::communicationsPhysAddAMR(Cellule **cellules, const int &dim, const int &lvl)
{
	Calcul_Parallele.communicationsVecteurAMR(cellules, "GPA", dim, lvl, m_numGPAGradC);
}

//***********************************************************************

int PAKCapillarite::getNumTransportAssocie() const
{
  return m_numTransportAssocie;
}

//***********************************************************************