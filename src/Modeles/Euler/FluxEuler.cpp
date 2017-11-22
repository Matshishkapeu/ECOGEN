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

#include <cmath>
#include "FluxEuler.h"

using namespace std;

FluxEuler fluxTempEuler;

//***********************************************************************

FluxEuler::FluxEuler(){}

//***********************************************************************

FluxEuler::~FluxEuler(){}

//***********************************************************************

void FluxEuler::afficheFlux() const
{
  cout << m_masse << " " << m_qdm.getX() << " " << m_energ << endl;
}

//***********************************************************************

void FluxEuler::ajoutFlux(double coefA, const int &nombrePhases)
{
    m_masse += coefA*fluxTempEuler.m_masse;
    m_qdm   += coefA*fluxTempEuler.m_qdm;
    m_energ += coefA*fluxTempEuler.m_energ;
}

//***********************************************************************

void FluxEuler::retireFlux(double coefA, const int &nombrePhases)
{
    m_masse -= coefA*fluxTempEuler.m_masse;
    m_qdm   -= coefA*fluxTempEuler.m_qdm;
    m_energ -= coefA*fluxTempEuler.m_energ;
}

//***********************************************************************

void FluxEuler::multiplie(double scalaire, const int &nombrePhases)
{
    m_masse *= scalaire;
    m_qdm   *= scalaire;
    m_energ *= scalaire;
}

//***********************************************************************

void FluxEuler::miseEnTampon(Cellule &cell, const int &nombrePhases)
{
  fluxTempEuler.construitCons(cell.getPhases(), nombrePhases, cell.getMelange());
}

//***********************************************************************

void FluxEuler::construitCons(Phase **phases, const int &nombrePhases, Melange *melange)
{
  Phase* phase(phases[0]);

  Eos *eos(phase->getEos());
  
  m_masse = phase->getDensite();
  m_qdm = m_masse*phase->getVitesse();
  m_energ = m_masse*phase->getEnergieTotale();
}

//***********************************************************************

void FluxEuler::construitPrim(Phase **phases, Melange *melange, const int &nombrePhases)
{
  double pression(0.), energieInterne(0.), energieTotale(0.), vitesseSon(0.);
  Phase* phase(phases[0]);
  Eos *eos(phase->getEos());

  phase->setDensite(m_masse);
  phase->setVitesse(m_qdm.getX() / m_masse, m_qdm.getY() / m_masse, m_qdm.getZ() / m_masse);
  //Mise a zero des petites vitesses
  if (fabs(phase->getU()) < 1.e-8) phase->setU(0.);
  if (fabs(phase->getV()) < 1.e-8) phase->setV(0.);
  if (fabs(phase->getW()) < 1.e-8) phase->setW(0.);

  energieTotale = m_energ / m_masse;
  phase->setEnergieTotale(energieTotale);
  energieInterne = energieTotale - 0.5*phase->getVitesse().normeCarre();
  phase->setEnergie(energieInterne);
  pression = eos->calculPression(m_masse, energieInterne);
  phase->setPression(pression);
  vitesseSon = eos->calculVitesseSon(m_masse, pression);
  phase->setVitesseSon(vitesseSon);
}

//***********************************************************************

void FluxEuler::miseAZero(const int &nombrePhases)
{
  m_masse = 0.;
  m_qdm   = 0.;
  m_energ = 0.;
}

//***********************************************************************

void FluxEuler::ajoutTuyere1D(Coord &normale, double const &surface, Cellule *cell, const int &nombrePhases)
{
  double coef = normale.getX()*surface / cell->getElement()->getVolume();
  Phase * phase;
  phase = cell->getPhase(0);

  m_qdm.setX(m_qdm.getX() - phase->getPression()*coef);
}

//***********************************************************************

void FluxEuler::retireTuyere1D(Coord &normale, double const &surface, Cellule *cell, const int &nombrePhases)
{
  double coef = normale.getX()*surface / cell->getElement()->getVolume();
  Phase * phase;
  phase = cell->getPhase(0);

  m_qdm.setX(m_qdm.getX() + phase->getPression()*coef);

}

//***********************************************************************

Coord FluxEuler::getQdm() const
{
  return m_qdm;
}

//***********************************************************************

double FluxEuler::getMasseMel() const
{
  return m_masse;
}

//***********************************************************************

double FluxEuler::getEnergieMel() const
{
  return m_energ;
}

//***********************************************************************

void FluxEuler::setCons(const Flux *cons, const int &nombrePhases)
{
  m_masse = cons->getMasseMel();
  m_qdm = cons->getQdm();
  m_energ = cons->getEnergieMel();
}

//***********************************************************************