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

#include "CondLimMurO2.h"

using namespace std;

//****************************************************************************

CondLimMurO2::CondLimMurO2() {}

//****************************************************************************

CondLimMurO2::CondLimMurO2(const CondLimMurO2& Source, const int lvl) : CondLimMur(Source, lvl)
{}

//****************************************************************************

CondLimMurO2::CondLimMurO2(int numPhysique) : CondLimMur(numPhysique)
{}

//****************************************************************************

CondLimMurO2::~CondLimMurO2()
{
  for (int k = 0; k < m_nombrePhases; k++) {
    delete m_vecPhasesPentes[k];
  }
  delete[] m_vecPhasesPentes;
  delete m_melangePentes;
  delete[] m_vecTransportsPentes;
}

//****************************************************************************

void CondLimMurO2::creeLimite(BordDeMaille **face)
{
  *face = new CondLimMurO2(*(this));
}

//***********************************************************************

void CondLimMurO2::allouePentes(const int &nombrePhases, const int &nombreTransports, int &allouePenteLocal)
{
  m_nombrePhases = nombrePhases;

  //Allocation des pentes des phases
  m_vecPhasesPentes = new Phase*[nombrePhases];
  //On attribut les phases a partir de la cellule a gauche (car cellule a droite inexistante pour les limites)
  //Necessaire car il faut connaitre le type de phase (ex: PhaseKapila, etc.))
  //Ensuite on met a zero toutes les pentes
  for (int k = 0; k < nombrePhases; k++) {
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
}

//***********************************************************************

void CondLimMurO2::calculPentes(const int &nombrePhases, const int &nombreTransports, Prim type)
{
  //Pentes des vitesses normales aux limites
  double distanceX, distanceY, distanceZ;
  distanceX = m_cellGauche->getElement()->distanceX(m_face);
  distanceY = m_cellGauche->getElement()->distanceY(m_face);
  distanceZ = m_cellGauche->getElement()->distanceZ(m_face);
  if (fabs(distanceX) > 1.e-8) {
    for (int k = 0; k < nombrePhases; k++) {
      m_vecPhasesPentes[k]->setU(m_cellGauche->getPhase(k, type)->getVitesse().getX() / distanceX);
    }
    m_melangePentes->setU(m_cellGauche->getMelange(type)->getVitesse().getX() / distanceX);
  }
  if (fabs(distanceY) > 1.e-8) {
    for (int k = 0; k < nombrePhases; k++) {
      m_vecPhasesPentes[k]->setV(m_cellGauche->getPhase(k, type)->getVitesse().getY() / distanceY);
    }
    m_melangePentes->setV(m_cellGauche->getMelange(type)->getVitesse().getY() / distanceY);
  }
  if (fabs(distanceZ) > 1.e-8) {
    for (int k = 0; k < nombrePhases; k++) {
      m_vecPhasesPentes[k]->setW(m_cellGauche->getPhase(k, type)->getVitesse().getZ() / distanceZ);
    }
    m_melangePentes->setW(m_cellGauche->getMelange(type)->getVitesse().getZ() / distanceZ);
  }
}

//***********************************************************************

void CondLimMurO2::resolRiemann(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type)
{
  cellGauche->copieVec(m_cellGauche->getPhases(type), m_cellGauche->getMelange(type), m_cellGauche->getTransports(type));

  //Calcul des distances bord de maille <-> cellules pour l extrapolation
  double distanceGauche(this->distance(m_cellGauche));
  //KS//FP// A voir avec Fabien comment faire ca bien pour qu'il n'y est aucun probleme en non-structure !
  //Probleme que la distance est une norme et donc on ne change pas le signe de la pente pour faire correctement l'extrapolation
  normale = m_face->getNormale();
  tangente = m_face->getTangente();
  binormale = m_face->getBinormale();
  if (normale.getX() < 0. || normale.getY() < 0. || normale.getZ() < 0.) { distanceGauche = -distanceGauche; }

  //Extrapolation gauche
  m_cellGauche->calculPentesLocalLimite(nombrePhases, nombreTransports, *this, limiteurGlobal, limiteurInterface);
  for (int k = 0; k < nombrePhases; k++) {
    cellGauche->getPhase(k)->extrapole(*pentesPhasesLocal1[k], distanceGauche);
  }
  cellGauche->getMelange()->extrapole(*pentesMelangeLocal1, distanceGauche);
  for (int k = 0; k < nombreTransports; k++) {
    cellGauche->getTransport(k).extrapole(pentesTransportLocal1[k], distanceGauche);
  }

  //Projection des vitesses sur repere attache a la face
  cellGauche->projection(normale, tangente, binormale, nombrePhases);
  //Calcul des variables etendus (Phases, Melange, PhysAdd)
  cellGauche->calculsEtendusPourRiemann(nombrePhases);

  //Probleme de Riemann
  double dxGauche(m_cellGauche->getElement()->getLCFL());
  dxGauche = dxGauche*pow(2., (double)m_lvl);
  this->resolRiemannLimite(*cellGauche, nombrePhases, dxGauche, dtMax);
  //Traitement des fonctions de transport (m_Sm connu : doit etre place apres l appel au Solveur de Riemann)
  if (nombreTransports > 0) { this->resolRiemannTransportLimite(*cellGauche, nombreTransports); }

  //Projection du flux sur le repere absolu
  m_mod->projectionRepereAbsolu(normale, tangente, binormale);
}

//***********************************************************************

Phase* CondLimMurO2::getPentesPhase(const int &numeroPhase) const
{
  return m_vecPhasesPentes[numeroPhase];
}

//***********************************************************************

Melange* CondLimMurO2::getPentesMelange() const
{
  return m_melangePentes;
}

//***********************************************************************

Transport* CondLimMurO2::getPentesTransport(const int &numeroTransport) const
{
  return &m_vecTransportsPentes[numeroTransport];
}


//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void CondLimMurO2::creerBordEnfant()
{
  m_bordsEnfants.push_back(new CondLimMurO2(*this, m_lvl + 1));
}

//****************************************************************************