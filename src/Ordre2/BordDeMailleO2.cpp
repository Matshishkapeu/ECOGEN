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

#include "BordDeMailleO2.h"
#include <iostream>

using namespace std;

Phase **pentesPhasesLocal1;
Phase **pentesPhasesLocal2;
Melange *pentesMelangeLocal1;
Melange *pentesMelangeLocal2;
double *pentesTransportLocal1;
double *pentesTransportLocal2;

//***********************************************************************

BordDeMailleO2::BordDeMailleO2() : BordDeMaille(), m_vecPhasesPentes(0), m_melangePentes(0), m_vecTransportsPentes(0)
 //m_BG1M(0), m_BG2M(0), m_BG3M(0), m_BG1P(0), m_BG2P(0), m_BG3P(0), m_BD1M(0),
 //m_BD2M(0), m_BD3M(0), m_BD1P(0), m_BD2P(0), m_BD3P(0),
 //m_betaG1M(1.), m_betaG2M(0.), m_betaG3M(0.), m_betaG1P(1.), m_betaG2P(0.), m_betaG3P(0.),
 //m_betaD1M(1.), m_betaD2M(0.), m_betaD3M(0.), m_betaD1P(1.), m_betaD2P(0.), m_betaD3P(0.),
 //m_distanceHGM(0.), m_distanceHGP(0.), m_distanceHDM(0.), m_distanceHDP(0.)
 {}

//***********************************************************************

BordDeMailleO2::BordDeMailleO2(int lvl) : BordDeMaille(lvl), m_vecPhasesPentes(0), m_melangePentes(0), m_vecTransportsPentes(0)
//m_BG1M(0), m_BG2M(0), m_BG3M(0), m_BG1P(0), m_BG2P(0), m_BG3P(0), m_BD1M(0),
//m_BD2M(0), m_BD3M(0), m_BD1P(0), m_BD2P(0), m_BD3P(0),
//m_betaG1M(1.), m_betaG2M(0.), m_betaG3M(0.), m_betaG1P(1.), m_betaG2P(0.), m_betaG3P(0.),
//m_betaD1M(1.), m_betaD2M(0.), m_betaD3M(0.), m_betaD1P(1.), m_betaD2P(0.), m_betaD3P(0.),
//m_distanceHGM(0.), m_distanceHGP(0.), m_distanceHDM(0.), m_distanceHDP(0.)
{}

//***********************************************************************

BordDeMailleO2::~BordDeMailleO2()
{
  for (int k = 0; k < m_nombrePhases; k++) {
    delete m_vecPhasesPentes[k];
  }
  delete[] m_vecPhasesPentes;
  delete m_melangePentes;
  delete[] m_vecTransportsPentes;
}

//***********************************************************************

void BordDeMailleO2::allouePentes(const int &nombrePhases, const int &nombreTransports, int &allouePenteLocal)
{
  m_nombrePhases = nombrePhases;

  //Allocation des pentes des phases
	m_vecPhasesPentes = new Phase*[nombrePhases];
  //On attribut les phases a partir de la cellule a gauche (car cellule a droite inexistante pour les limites)
  //Necessaire car il faut connaitre le type de phase (ex: PhaseKapila, etc.))
  //Ensuite on met a zero toutes les pentes
  for(int k = 0; k < nombrePhases; k++){
    m_cellGauche->getPhase(k)->alloueEtCopiePhase(&m_vecPhasesPentes[k]);
		m_vecPhasesPentes[k]->miseAZero();
  }
  m_cellGauche->getMelange()->alloueEtCopieMelange(&m_melangePentes);
  m_melangePentes->miseAZero();

	//Allocation des pentes sur transports
	m_vecTransportsPentes = new Transport[nombreTransports];
	for (int k = 0; k < nombreTransports; k++) {
		m_vecTransportsPentes[k].setValeur(0.);
	}

  //Allocation des variables externes
  if (allouePenteLocal < 1) {
    pentesPhasesLocal1 = new Phase*[nombrePhases];
    pentesPhasesLocal2 = new Phase*[nombrePhases];
    for (int k = 0; k < nombrePhases; k++) {
      m_cellGauche->getPhase(k)->alloueEtCopiePhase(&pentesPhasesLocal1[k]);
      m_cellGauche->getPhase(k)->alloueEtCopiePhase(&pentesPhasesLocal2[k]);
      pentesPhasesLocal1[k]->miseAZero();
      pentesPhasesLocal2[k]->miseAZero();
    }

    m_cellGauche->getMelange()->alloueEtCopieMelange(&pentesMelangeLocal1);
    m_cellGauche->getMelange()->alloueEtCopieMelange(&pentesMelangeLocal2);
    pentesMelangeLocal1->miseAZero();
    pentesMelangeLocal2->miseAZero();

		pentesTransportLocal1 = new double[nombreTransports];
		pentesTransportLocal2 = new double[nombreTransports];
		for (int k = 0; k < nombreTransports; k++) {
			pentesTransportLocal1[k] = 0.;
			pentesTransportLocal2[k] = 0.;
		}

    allouePenteLocal = 1;
  }
}

//***********************************************************************

void BordDeMailleO2::calculPentes(const int &nombrePhases, const int &nombreTransports, Prim type)
{
  if (m_bordsEnfants.size() == 0) {
    //Distance entre les deux mailles en contact
    double distance(m_cellGauche->distance(m_cellDroite));
    //Attribution gauche/droite
    //Si la cellule gauche ou droite est de niveau inferieur a "lvl", on ne prend pas "type" mais vecPhases (ca evite de prendre vecPhaseO2 alors qu'on ne l'a pas).
    Prim typeGauche = type;
    Prim typeDroite = type;
    if (m_cellGauche->getLvl() < m_lvl) { typeGauche = vecPhases; }
    if (m_cellDroite->getLvl() < m_lvl) { typeDroite = vecPhases; }

    for (int k = 0; k < nombrePhases; k++) {
			m_vecPhasesPentes[k]->calculPentesPhase(*m_cellGauche->getPhase(k, typeGauche), *m_cellDroite->getPhase(k, typeDroite), distance);
    }
    m_melangePentes->calculPentesMelange(*m_cellGauche->getMelange(typeGauche), *m_cellDroite->getMelange(typeDroite), distance);
    for (int k = 0; k < nombreTransports; k++) {
      m_vecTransportsPentes[k].calculPentesTransport(m_cellGauche->getTransport(k, typeGauche).getValeur(), m_cellDroite->getTransport(k, typeDroite).getValeur(), distance);
    }
  }
}

//***********************************************************************

void BordDeMailleO2::calculFlux(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type)
{
  // Quand on fait le premier calculFlux (donc avec vecPhases) on n'incremente pas m_cons pour les mailles de niveau different (inferieur) de "lvl".
  // Sinon va veut dire qu on l ajoute pour les 2 calculFlux sans le remettre a zero entre les deux, donc 2 fois plus de flux que ce que l on veut.
  this->resolRiemann(nombrePhases, nombreTransports, dtMax, limiteurGlobal, limiteurInterface, type);

  switch (type) {
  case vecPhases:
    if (m_cellGauche->getLvl() == m_cellDroite->getLvl()) {     //CoefAMR = 1 pour les deux
      this->ajoutFlux(nombrePhases, nombreTransports, 1.);      //Ajout du flux sur maille droite
      this->retireFlux(nombrePhases, nombreTransports, 1.);     //Retrait du flux sur maille gauche
    }
    else if (m_cellGauche->getLvl() > m_cellDroite->getLvl()) { //CoefAMR = 1 pour la gauche et on n'ajoute rien sur la maille droite
      this->retireFlux(nombrePhases, nombreTransports, 1.);     //Retrait du flux sur maille gauche
    }
    else {                                                      //CoefAMR = 1 pour la droite et on ne retire rien sur la maille gauche
      this->ajoutFlux(nombrePhases, nombreTransports, 1.);      //Ajout du flux sur maille droite
    }
    break;

  case vecPhasesO2:
    if (m_cellGauche->getLvl() == m_cellDroite->getLvl()) {     //CoefAMR = 1 pour les deux
      this->ajoutFlux(nombrePhases, nombreTransports, 1.);      //Ajout du flux sur maille droite
      this->retireFlux(nombrePhases, nombreTransports, 1.);     //Retrait du flux sur maille gauche
    }
    else if (m_cellGauche->getLvl() > m_cellDroite->getLvl()) { //CoefAMR = 1 pour la gauche et 0.5 pour la droite
      this->ajoutFlux(nombrePhases, nombreTransports, 0.5);     //Ajout du flux sur maille droite
      this->retireFlux(nombrePhases, nombreTransports, 1.);     //Retrait du flux sur maille gauche
    }
    else {                                                      //CoefAMR = 1 pour la droite et 0.5 pour la gauche
      this->ajoutFlux(nombrePhases, nombreTransports, 1.);      //Ajout du flux sur maille droite
      this->retireFlux(nombrePhases, nombreTransports, 0.5);    //Retrait du flux sur maille gauche
    }
    break;

  default: break;
  }
}

//***********************************************************************

void BordDeMailleO2::resolRiemann(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type)
{
  //Si la cellule gauche ou droite est de niveau inferieur a "lvl", on ne prend pas "type" mais vecPhases (ca evite de prendre vecPhaseO2 alors qu'on ne l'a pas).
  if (m_cellGauche->getLvl() == m_lvl) { cellGauche->copieVec(m_cellGauche->getPhases(type), m_cellGauche->getMelange(type), m_cellGauche->getTransports(type)); }
  else { cellGauche->copieVec(m_cellGauche->getPhases(vecPhases), m_cellGauche->getMelange(vecPhases), m_cellGauche->getTransports(vecPhases)); }
  
  if (m_cellDroite->getLvl() == m_lvl) { cellDroite->copieVec(m_cellDroite->getPhases(type), m_cellDroite->getMelange(type), m_cellDroite->getTransports(type)); }
  else { cellDroite->copieVec(m_cellDroite->getPhases(vecPhases), m_cellDroite->getMelange(vecPhases), m_cellDroite->getTransports(vecPhases)); }

  //Calcul des distances bord de maille <-> cellules pour l extrapolation
  double distanceGauche(this->distance(m_cellGauche));
  double distanceDroite(this->distance(m_cellDroite));

  //Extrapolation gauche
  m_cellGauche->calculPentesLocal(nombrePhases, nombreTransports, *this, limiteurGlobal, limiteurInterface);
  for (int k = 0; k < nombrePhases; k++) {
    cellGauche->getPhase(k)->extrapole(*pentesPhasesLocal1[k], distanceGauche);
  }
  cellGauche->getMelange()->extrapole(*pentesMelangeLocal1, distanceGauche);
	for (int k = 0; k < nombreTransports; k++) {
		cellGauche->getTransport(k).extrapole(pentesTransportLocal1[k], distanceGauche);
	}

  //Extrapolation droite
  m_cellDroite->calculPentesLocal(nombrePhases, nombreTransports, *this, limiteurGlobal, limiteurInterface);
  for (int k = 0; k < nombrePhases; k++) {
    pentesPhasesLocal1[k]->changeSigne(); //On doit soustraire les pentes a droite
    cellDroite->getPhase(k)->extrapole(*pentesPhasesLocal1[k], distanceDroite);
  }
  pentesMelangeLocal1->changeSigne();
  cellDroite->getMelange()->extrapole(*pentesMelangeLocal1, distanceDroite);
	for (int k = 0; k < nombreTransports; k++) {
		pentesTransportLocal1[k] = -pentesTransportLocal1[k];
		cellDroite->getTransport(k).extrapole(pentesTransportLocal1[k], distanceDroite);
	}

  //Projection des vitesses sur repere attache a la face
  normale = m_face->getNormale();
  tangente = m_face->getTangente();
  binormale = m_face->getBinormale();
  cellGauche->projection(normale, tangente, binormale, nombrePhases);
  cellDroite->projection(normale, tangente, binormale, nombrePhases);

  //Calcul des variables etendus (Phases, Melange, PhysAdd)
  cellGauche->calculsEtendusPourRiemann(nombrePhases);
  cellDroite->calculsEtendusPourRiemann(nombrePhases);

  //Probleme de Riemann
  double dxGauche(m_cellGauche->getElement()->getLCFL());
  double dxDroite(m_cellDroite->getElement()->getLCFL());
  dxGauche = dxGauche*pow(2., (double)m_lvl);
  dxDroite = dxDroite*pow(2., (double)m_lvl);
  m_mod->resolRiemannInterne(*cellGauche, *cellDroite, nombrePhases, dxGauche, dxDroite, dtMax);
  //Traitement des fonctions de transport (m_Sm connu : doit etre place apres l appel au Solveur de Riemann)
  if (nombreTransports > 0) { m_mod->resolRiemannTransportInterne(*cellGauche, *cellDroite, nombreTransports); }

  //Projection du flux sur le repere absolu
  m_mod->projectionRepereAbsolu(normale, tangente, binormale);
}

//***********************************************************************

Phase* BordDeMailleO2::getPentesPhase(const int &numeroPhase) const
{
  return m_vecPhasesPentes[numeroPhase];
}

//***********************************************************************

Melange* BordDeMailleO2::getPentesMelange() const
{
  return m_melangePentes;
}

//***********************************************************************

Transport* BordDeMailleO2::getPentesTransport(const int &numeroTransport) const
{
  return &m_vecTransportsPentes[numeroTransport];
}

//***********************************************************************

//Cellule * BordDeMailleO2::getB(BO2 B) const 
//{
//  switch (B){
//  case BG1M: return m_BG1M; break;
//  case BG2M: return m_BG2M; break;
//  case BG1P: return m_BG1P; break;
//  case BG2P: return m_BG2P; break;
//  case BD1M: return m_BD1M; break;
//  case BD2M: return m_BD2M; break;
//  case BD1P: return m_BD1P; break;
//  case BD2P: return m_BD2P; break;
//  default: Erreurs::messageErreur("probleme enum non connu dans BorDeMailleO2::getB"); return 0; break;
//  }
//}

//***********************************************************************

//double BordDeMailleO2::getBeta(betaO2 beta) const 
//{
//  switch (beta) {
//  case betaG1M: return m_betaG1M; break;
//  case betaG2M: return m_betaG2M; break;
//  case betaG1P: return m_betaG1P; break;
//  case betaG2P: return m_betaG2P; break;
//  case betaD1M: return m_betaD1M; break;
//  case betaD2M: return m_betaD2M; break;
//  case betaD1P: return m_betaD1P; break;
//  case betaD2P: return m_betaD2P; break;
//  default: Erreurs::messageErreur("probleme enum non connu dans BorDeMailleO2::getBeta"); return 0; break;
//  }
//}

//***********************************************************************

//double BordDeMailleO2::getDistanceH(distanceHO2 dist) const
//{
//  switch (dist) {
//  case distanceHGM: return m_distanceHGM; break;
//  case distanceHGP: return m_distanceHGP; break;
//  case distanceHDM: return m_distanceHDM; break;
//  case distanceHDP: return m_distanceHDP; break;
//  default: Erreurs::messageErreur("probleme enum non connu dans BorDeMailleO2::getDistanceH"); return 0; break;
//  }
//}

//***********************************************************************

//void BordDeMailleO2::setB(BO2 B, Cellule *cellule)
//{
//  switch (B) {
//  case BG1M: m_BG1M = cellule; break;
//  case BG2M: m_BG2M = cellule; break;
//  case BG1P: m_BG1P = cellule; break;
//  case BG2P: m_BG2P = cellule; break;
//  case BD1M: m_BD1M = cellule; break;
//  case BD2M: m_BD2M = cellule; break;
//  case BD1P: m_BD1P = cellule; break;
//  case BD2P: m_BD2P = cellule; break;
//  default: Erreurs::messageErreur("probleme enum non connu dans BorDeMailleO2::setB"); break;
//  }
//}

//***********************************************************************

//void BordDeMailleO2::setBeta(betaO2 beta, double &valeur)
//{
//  switch (beta) {
//  case betaG1M: m_betaG1M = valeur; break;
//  case betaG2M: m_betaG2M = valeur; break;
//  case betaG1P: m_betaG1P = valeur; break;
//  case betaG2P: m_betaG2P = valeur; break;
//  case betaD1M: m_betaD1M = valeur; break;
//  case betaD2M: m_betaD2M = valeur; break;
//  case betaD1P: m_betaD1P = valeur; break;
//  case betaD2P: m_betaD2P = valeur; break;
//  default: Erreurs::messageErreur("probleme enum non connu dans BorDeMailleO2::setB"); break;
//  }
//}

//***********************************************************************

//void BordDeMailleO2::setDistanceH(distanceHO2 dist, double &valeur) 
//{
//  switch (dist) {
//  case distanceHGM: m_distanceHGM = valeur; break;
//  case distanceHGP: m_distanceHGP = valeur; break;
//  case distanceHDM: m_distanceHDM = valeur; break;
//  case distanceHDP: m_distanceHDP = valeur; break;
//  default: Erreurs::messageErreur("probleme enum non connu dans BorDeMailleO2::setDistanceH"); break;
//  }
//}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void BordDeMailleO2::creerBordEnfant()
{
  m_bordsEnfants.push_back(new BordDeMailleO2(m_lvl + 1));
}

//***********************************************************************

void BordDeMailleO2::creerBordEnfantInterne(const int &lvl, vector<BordDeMaille*> *bordsEnfantsInternes)
{
  (*bordsEnfantsInternes).push_back(new BordDeMailleO2(lvl + 1));
}

//***********************************************************************