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
#include "FluxEulerHomogene.h"

using namespace std;

FluxEulerHomogene fluxTempEulerHomogene;

//***********************************************************************

FluxEulerHomogene::FluxEulerHomogene() {}

//***********************************************************************

FluxEulerHomogene::FluxEulerHomogene(ModEulerHomogene *modele) : m_modele(modele)
{}

//***********************************************************************

FluxEulerHomogene::~FluxEulerHomogene(){}

//***********************************************************************

void FluxEulerHomogene::afficheFlux() const
{
  cout << m_masse << " " << m_qdm.getX() << " " << m_energ << endl;
}

//***********************************************************************

void FluxEulerHomogene::ajoutFlux(double coefA, const int &nombrePhases)
{
    m_masse += coefA*fluxTempEulerHomogene.m_masse;
    m_qdm   += coefA*fluxTempEulerHomogene.m_qdm;
    m_energ += coefA*fluxTempEulerHomogene.m_energ;
}

//***********************************************************************

void FluxEulerHomogene::retireFlux(double coefA, const int &nombrePhases)
{
    m_masse -= coefA*fluxTempEulerHomogene.m_masse;
    m_qdm   -= coefA*fluxTempEulerHomogene.m_qdm;
    m_energ -= coefA*fluxTempEulerHomogene.m_energ;
}

//***********************************************************************

void FluxEulerHomogene::multiplie(double scalaire, const int &nombrePhases)
{
    m_masse *= scalaire;
    m_qdm   *= scalaire;
    m_energ *= scalaire;
}

//***********************************************************************

void FluxEulerHomogene::miseEnTampon(Cellule &cell, const int &nombrePhases)
{
  fluxTempEulerHomogene.construitCons(cell.getPhases(), nombrePhases, cell.getMelange());
}

//***********************************************************************

void FluxEulerHomogene::construitCons(Phase **phases, const int &nombrePhases, Melange *melange)
{
  double energieInterne(0.);
  double rhok, alphak, ek;

  Phase *phase(0);
  m_masse = 0.;
  m_energ = 0.;

  for (int k = 0; k < nombrePhases; k++)
  {
    phase = phases[k];
    //Calcul masse volumique de melange
    alphak = phase->getAlpha();
    rhok = phase->getDensite();
    m_masse += alphak * rhok;
    //Calcul energie interne massique de melange
    ek = phase->getEos()->calculEnergie(rhok, phase->getPression());
    energieInterne += alphak*rhok*ek;
  }
  m_qdm = m_masse*melange->getVitesse();
  m_energ = energieInterne + 0.5*m_masse*melange->getVitesse().normeCarre();
}

//***********************************************************************

void FluxEulerHomogene::construitPrim(Phase **phases, Melange *melange, const int &nombrePhases)
{
  double pression, energieInterne, energieTotale(0.), vitesseSon;
  Phase* phase(0);
  
  int liq(m_modele->m_liq), vap(m_modele->m_vap);
  double rhoe(m_energ - 0.5*m_qdm.normeCarre() / m_masse); //Q// Pourquoi diviser par rho si c est rhoe ???...

  //Recherche de la pression par procede iteratif sur conservation de l'energie (e=Somme(Ykek))
  int iteration(0);
  pression = phases[1]->getPression(); //estimation pression precedente
  double f(0.), df(1.);
  double alphaVap, dalphaVap, rhoLiq, drhoLiq, rhoVap, drhoVap, rhoeLiq, drhoeLiq, rhoeVap, drhoeVap, Tsat, dTsat;
  do {
    pression -= f / df; iteration++;
    if (iteration > 50) {
      erreurs.push_back(Erreurs("nombre iterations trop grand dans construitPrim EulerHomogene", __FILE__, __LINE__));
      break;
    }
    Tsat = melange->calculTsat(phases[liq]->getEos(), phases[vap]->getEos(), pression, &dTsat);
    rhoVap = phases[vap]->getEos()->calculDensiteSaturation(pression, Tsat, dTsat, &drhoVap);
    rhoLiq = phases[liq]->getEos()->calculDensiteSaturation(pression, Tsat, dTsat, &drhoLiq);
    alphaVap = (m_masse - rhoLiq) / (rhoVap - rhoLiq);
    dalphaVap = (-drhoLiq*(rhoVap - rhoLiq) - (m_masse - rhoLiq)*(drhoVap - drhoLiq)) / ((rhoVap - rhoLiq)*(rhoVap - rhoLiq));
    rhoeVap = phases[vap]->getEos()->calculRhoEnergieSaturation(pression, rhoVap, drhoVap, &drhoeVap);
    rhoeLiq = phases[liq]->getEos()->calculRhoEnergieSaturation(pression, rhoLiq, drhoLiq, &drhoeLiq);
    f = rhoe - alphaVap*(rhoeVap - rhoeLiq) - rhoeLiq;
    df = -dalphaVap*(rhoeVap - rhoeLiq) - alphaVap*(drhoeVap - drhoeLiq) - drhoeLiq;
  } while (fabs(f/rhoe) > 1e-10);

//  cout << rhoVap << " " << rhoLiq << endl;
//  cout << pression << " " << Tsat << endl;
//  exit(0);

  //Variables des phases
  phases[liq]->setAlpha(1. - alphaVap);
  phases[liq]->setDensite(rhoLiq);

  phases[vap]->setAlpha(alphaVap);
  phases[vap]->setDensite(rhoVap);
  
  energieTotale = m_energ / m_masse;
  melange->setVitesse(m_qdm / m_masse);
  //Mise a zero des petites vitesses
  if (fabs(melange->getU()) < 1.e-8) melange->setU(0.);
  if (fabs(melange->getV()) < 1.e-8) melange->setV(0.);
  if (fabs(melange->getW()) < 1.e-8) melange->setW(0.);
  for (int k = 0; k < 2; k++) {
    phases[k]->setPression(pression);
    phases[k]->setEnergieTotale(energieTotale);
    energieInterne = energieTotale - 0.5*melange->getVitesse().normeCarre();
    phases[k]->setEnergie(energieInterne);
    vitesseSon = phases[k]->getEos()->calculVitesseSon(m_masse, pression);
    phases[k]->setVitesseSon(vitesseSon);
  }

  //Variables de melange
  melange->calculGrandeursMelange(phases, nombrePhases);
  melange->setEnergieTotale(energieTotale);
}

//***********************************************************************

void FluxEulerHomogene::miseAZero(const int &nombrePhases)
{
  m_masse = 0.;
  m_qdm   = 0.;
  m_energ = 0.;
}

//***********************************************************************

void FluxEulerHomogene::ajoutTuyere1D(Coord &normale, double const &surface, Cellule *cell, const int &nombrePhases)
{
  double coef = normale.getX()*surface / cell->getElement()->getVolume();
  Phase * phase;
  phase = cell->getPhase(0);

  m_qdm.setX(m_qdm.getX() - phase->getPression()*coef);
}
//***********************************************************************

void FluxEulerHomogene::retireTuyere1D(Coord &normale, double const &surface, Cellule *cell, const int &nombrePhases)
{
  double coef = normale.getX()*surface / cell->getElement()->getVolume();
  Phase * phase;
  phase = cell->getPhase(0);

  m_qdm.setX(m_qdm.getX() + phase->getPression()*coef);

}

//***********************************************************************

Coord FluxEulerHomogene::getQdm() const
{
  return m_qdm;
}

//***********************************************************************

double FluxEulerHomogene::getMasseMel() const
{
  return m_masse;
}

//***********************************************************************

double FluxEulerHomogene::getEnergieMel() const
{
  return m_energ;
}

//***********************************************************************

void FluxEulerHomogene::setCons(const Flux *cons, const int &nombrePhases)
{
  m_masse = cons->getMasseMel();
  m_qdm = cons->getQdm();
  m_energ = cons->getEnergieMel();
}

//***********************************************************************