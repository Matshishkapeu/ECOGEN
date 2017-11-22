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
#include <algorithm>
#include "FluxKapila.h"
#include "../Melange.h"

using namespace std;

FluxKapila *fluxTempKapila;
FluxKapila *sourceConsKap;

//***********************************************************************

FluxKapila::FluxKapila(){}

//***********************************************************************

FluxKapila::FluxKapila(const int &nombrePhases)
{
  m_alpha = new double[nombrePhases];
  m_masse = new double[nombrePhases];
  m_energ = new double[nombrePhases];
}

//***********************************************************************

FluxKapila::~FluxKapila()
{
  delete[] m_alpha;
  delete[] m_masse;
  delete[] m_energ;
}

//***********************************************************************

void FluxKapila::afficheFlux() const
{
  cout << m_masse << " " << m_qdm.getX() << " " << m_energ << endl;
}

//***********************************************************************

void FluxKapila::ajoutFlux(double coefA, const int &nombrePhases)
{
  for (int k = 0; k < nombrePhases; k++) {
    m_alpha[k] += coefA*fluxTempKapila->m_alpha[k];
    m_masse[k] += coefA*fluxTempKapila->m_masse[k];
    m_energ[k] += coefA*fluxTempKapila->m_energ[k];
  }
  m_qdm += coefA*fluxTempKapila->m_qdm;
  m_energMelange += coefA*fluxTempKapila->m_energMelange;
}

//***********************************************************************

void FluxKapila::retireFlux(double coefA, const int &nombrePhases)
{
  for (int k = 0; k < nombrePhases; k++) {
    m_alpha[k] -= coefA*fluxTempKapila->m_alpha[k];
    m_masse[k] -= coefA*fluxTempKapila->m_masse[k];
    m_energ[k] -= coefA*fluxTempKapila->m_energ[k];
  }
  m_qdm -= coefA*fluxTempKapila->m_qdm;
  m_energMelange -= coefA*fluxTempKapila->m_energMelange;
}

//***********************************************************************

void FluxKapila::multiplie(double scalaire, const int &nombrePhases)
{
    for(int k=0;k<nombrePhases;k++)
    {
      m_alpha[k] *= scalaire;
      m_masse[k] *= scalaire;
      m_energ[k] *= scalaire;
    }
    m_qdm  *= scalaire;
    m_energMelange *= scalaire;
}

//***********************************************************************

void FluxKapila::miseEnTampon(Cellule &cell, const int &nombrePhases)
{
  fluxTempKapila->construitCons(cell.getPhases(), nombrePhases, cell.getMelange());
}

//***********************************************************************

void FluxKapila::construitCons(Phase **phases, const int &nombrePhases, Melange *melange)
{
	double energieInterne(0.);
	Phase *phase(0);
  for (int k = 0; k < nombrePhases; k++)
	{
		phase = phases[k];
		BO->ak[k] = phase->getAlpha();
		BO->rhok[k] = phase->getDensite();
		m_alpha[k] = BO->ak[k];
		m_masse[k] = BO->ak[k] * BO->rhok[k];
		//Calcul energie totale phase k
		energieInterne = BO->eos[k]->calculEnergie(BO->rhok[k], phase->getPression());
		m_energ[k] = BO->ak[k] * BO->rhok[k] * energieInterne;
	}
	m_qdm = melange->getDensite()*melange->getVitesse();
  m_energMelange = melange->getDensite()*melange->getEnergieTotale();
}

//***********************************************************************

void FluxKapila::construitPrim(Phase **phases, Melange *melange, const int &nombrePhases)
{
  double pression(0.), energieInterne(0.), rhoMel(0.);

  //Verification amont et correctif sur alpha et masse si necessaire (maintien de l ordre 2)
  double un(0.);
  for (int k = 0; k < nombrePhases; k++) {
    if (m_alpha[k] <= 1.e-10) m_alpha[k] = 1e-9;
    if (m_alpha[k] >= 1.-1.e-10) m_alpha[k] = 1.0 - 1e-9;
    un += m_alpha[k];
  }
  for (int k = 0; k < nombrePhases; k++) { m_alpha[k] /= un; }

  //Variables de chaque phase et melange
  Phase *phase(0);
  for (int k = 0; k < nombrePhases; k++) {
      phase = phases[k];
      rhoMel = rhoMel + m_masse[k];
      phase->setAlpha(m_alpha[k]);
      phase->setDensite(m_masse[k]/m_alpha[k]);
      //Calcul Pression
      energieInterne = m_energ[k]/m_masse[k];
      pression = BO->eos[k]->calculPression(phase->getDensite(),energieInterne);
      phase->setPression(pression);
  }
  melange->setVitesse(m_qdm.getX() / rhoMel, m_qdm.getY() / rhoMel, m_qdm.getZ() / rhoMel);
  //Mise a zero des petites vitesses
  if (fabs(melange->getU()) < 1.e-8) melange->setU(0.);
  if (fabs(melange->getV()) < 1.e-8) melange->setV(0.);
  if (fabs(melange->getW()) < 1.e-8) melange->setW(0.);
  for (int k = 0; k < nombrePhases; k++) {
    phases[k]->calculsEtendusPhases(melange->getVitesse());
  }
  melange->calculGrandeursMelange(phases, nombrePhases);
  //On reconstruit l energie totale par la valeur issue de l equation d energie totale.
  double energieTotaleMelange(0.);
  energieTotaleMelange = m_energMelange / rhoMel;
  melange->setEnergieTotale(energieTotaleMelange);
}

//***********************************************************************

void FluxKapila::miseAZero(const int &nombrePhases)
{
  for(int k=0;k<nombrePhases;k++){
    m_alpha[k] = 0.;
    m_masse[k] = 0.;
    m_energ[k] = 0.;
  }
  m_qdm = 0.;
  m_energMelange = 0.;
}

//***********************************************************************

void FluxKapila::miseAZeroFluxTemp(const int &nombrePhases)
{
  for (int k = 0; k<nombrePhases; k++) {
    fluxTempKapila->m_alpha[k] = 0.;
    fluxTempKapila->m_masse[k] = 0.;
    fluxTempKapila->m_energ[k] = 0.;
  }
  fluxTempKapila->m_qdm = 0.;
  fluxTempKapila->m_energMelange = 0.;
}

//***********************************************************************

void FluxKapila::ajoutNonCons(double coefA, const Cellule *cell, const int &nombrePhases)
{
  Phase *phase;
  for(int k=0;k<nombrePhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] += -coefA*phase->getAlpha()*fluxTempKapila->m_sM;
    m_energ[k] += coefA*phase->getAlpha()*phase->getPression()*fluxTempKapila->m_sM;
  }
}

//***********************************************************************

void FluxKapila::retireNonCons(double coefA, const Cellule *cell, const int &nombrePhases)
{
  Phase *phase;
  for(int k=0;k<nombrePhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] -= -coefA*phase->getAlpha()*fluxTempKapila->m_sM;
    m_energ[k] -= coefA*phase->getAlpha()*phase->getPression()*fluxTempKapila->m_sM;
  }
}

//***********************************************************************

void FluxKapila::relaxPressions(Cellule *cell, const int &nombrePhases, Prim type) const
{
  double drho(0.), dalpha(0.);
  Phase *phase(0);
  //Recuperation grandeurs initiales
  double pStar(0.);
  for (int k = 0; k < nombrePhases; k++)
  {
    phase = cell->getPhase(k, type);
    BO->ak[k] = phase->getAlpha();
    BO->pk[k] = phase->getPression();
    BO->rhok[k] = phase->getDensite();
    pStar += BO->ak[k] * BO->pk[k];
    phase->verifiePhase();
  }
  
  //Boucle iterative pour determination de la pression relaxee
  int iteration(0);
  double f(0.), df(1.);
  do {
    pStar -= f / df; iteration++;
    if (iteration > 50) {
      erreurs.push_back(Erreurs("nombre iterations trop grand dans relaxation Pression", __FILE__, __LINE__));
      break;
      //Erreurs::messageErreur("nombre iterations trop grand dans relaxation Pression");
    }
    //Verification pression physique ?
    for (int k = 0; k < nombrePhases; k++) { BO->eos[k]->verifieEtCorrigePression(pStar); }
    f = -1.; df = 0.;
    for (int k = 0; k<nombrePhases; k++)
    {
      BO->rhokS[k] = BO->eos[k]->calculDensiteIsentrope(BO->pk[k], BO->rhok[k], pStar, &drho);
      BO->akS[k] = BO->ak[k] * BO->rhok[k] / BO->rhokS[k];
      dalpha = BO->ak[k] * BO->rhok[k] * drho / (BO->rhokS[k] * BO->rhokS[k]);
      f += BO->akS[k];
      df -= dalpha;
    }
  } while (fabs(f)>1e-10);

  //Modification de la cellule par valeur relaxee en pression
  for (int k = 0; k<nombrePhases; k++)
  {
    phase = cell->getPhase(k, type);
    phase->setAlpha(BO->akS[k]);
    phase->setDensite(BO->rhokS[k]);
    phase->setPression(pStar);
  }
  cell->getMelange(type)->setPression(pStar);
}

//***********************************************************************

void FluxKapila::relaxPTMu(Cellule *cell, const int &nombrePhases, Prim type) const
{
  Phase *phase(0);

  int liq(1), vap(0); //couple liq-vap temporaire phase 0 et 1

  if (nombrePhases > 2) erreurs.push_back(Erreurs("evaporation a + de deux phases non prevue dans relaxPTMu", __FILE__, __LINE__));

  //Recuperation grandeurs initiales
  double pStar(0.), Tsat;
  for (int k = 0; k<nombrePhases; k++)
  {
    phase = cell->getPhase(k, type);
    BO->ak[k] = phase->getAlpha();
    BO->pk[k] = phase->getPression();
    BO->rhok[k] = phase->getDensite();
    pStar += BO->ak[k] * BO->pk[k];
    phase->verifiePhase();
  }
  //cell->calculsEtendus(nombrePhases);
  double rho = cell->getMelange()->getDensite();
  double rhoe = rho*cell->getMelange()->getEnergie();

  //Determination de la Temperature de saturation
  double dTsat(0.);
  Tsat = calculTsat(pStar, nombrePhases, &dTsat);
  
  //Pas evap ?
  double TL = BO->eos[liq]->calculTemperature(BO->rhok[liq],pStar);
  double TV = BO->eos[vap]->calculTemperature(BO->rhok[vap], pStar);
  if (TL < Tsat) return;

  //Boucle iterative pour determination de la pression relaxee
  double rhoLSat, rhoVSat, drhoLSat, drhoVSat;
  double aLSat, aVSat, daLSat, daVSat;
  double rhoeLSat, rhoeVSat, drhoeLSat, drhoeVSat;
  int iteration(0);
  double f(0.), df(1.);
  do {
    pStar -= f / df; iteration++;
    if (iteration > 50) {
      erreurs.push_back(Erreurs("nombre iterations trop grand dans relaxPTMu", __FILE__, __LINE__));
      cout << "info cellule prolbematique" << endl;
      cout << "Liq " << TL << " " << BO->rhok[liq] << " " << cell->getMelange()->getPression() << endl;
      cout << "Vap " << TV << " " << BO->rhok[vap] << " " << cell->getMelange()->getPression() << endl;
      cout << Tsat << " " << pStar << endl;
      break;
    }
    //Verification pression physique ?
    for (int k = 0; k < nombrePhases; k++) { BO->eos[k]->verifieEtCorrigePression(pStar); }
    //Determination des masses volumiques du couple liquide-vapeur
    Tsat = calculTsat(pStar, nombrePhases, &dTsat);
    rhoLSat = BO->eos[liq]->calculDensiteSaturation(pStar, Tsat, dTsat, &drhoLSat);
    rhoVSat = BO->eos[vap]->calculDensiteSaturation(pStar, Tsat, dTsat, &drhoVSat);
    //Limites pour les densites
    if (rhoLSat <= rho) {
      rhoLSat = rho + 1e-6;
      drhoLSat = 0.;
    }
    if (rhoVSat >= rho) {
      rhoVSat = rho - 1e-6;
      drhoVSat = 0.;
    }

    //Determination des fraction volumique couple liquide-vapeur
    aLSat = (rho - rhoVSat) / (rhoLSat - rhoVSat);
    daLSat = (- drhoVSat*(rhoLSat - rhoVSat)- (rho - rhoVSat)*(drhoLSat - drhoVSat)) / ((rhoLSat - rhoVSat)*(rhoLSat - rhoVSat));
    aVSat = 1. - aLSat; //KS//FP// Sert a rien celui là (mais laisser pour l'instant si objectif a changer)
    daVSat = -daLSat;
    //Limites pour les fractions volumiques
    if (aLSat <= 0.) aLSat = 1e-8;
    if (aLSat >= 1.) aLSat = 1.-1e-8;
    aVSat = 1. - aLSat;

    f = rhoe; df = 0.;
    //Couple liquide vapeur uniquement
    rhoeLSat = BO->eos[liq]->calculRhoEnergieSaturation(pStar, rhoLSat, drhoLSat, &drhoeLSat);
    rhoeVSat = BO->eos[vap]->calculRhoEnergieSaturation(pStar, rhoVSat, drhoVSat, &drhoeVSat);
    f -= (aLSat*rhoeLSat + aVSat*rhoeVSat);
    df -= (daLSat*rhoeLSat + aLSat*drhoeLSat + daVSat*rhoeVSat + aVSat*drhoeVSat);
    //cout << iteration << " " << pStar << " " << rho << " " << rhoLSat << " " << rhoVSat << " " << endl;
    //cout <<iteration<< " "<< pStar << " " << Tsat << " " << f << " "<< aLSat<< " "<< aVSat << endl;
    f /= rhoe;
    df /= rhoe;
  } while (fabs(f)>1e-10);

  //Mise a jour grandeurs des phases
  phase = cell->getPhase(liq);
  phase->setAlpha(aLSat);
  phase->setDensite(rhoLSat);
  phase->setPression(pStar);

  phase = cell->getPhase(vap);
  phase->setAlpha(aVSat);
  phase->setDensite(rhoVSat);
  phase->setPression(pStar);

  //cell->calculsEtendus(nombrePhases);
}

//***********************************************************************

double FluxKapila::calculTsat(const double &pression, const int &nombrePhases, double *dTsat) const
{
  //FP//TODO// ajouter un garde fou si les loi d etat ne sont pas adaptees
  int liq(1), vap(0); //couple liq-vap temporaire phase 0 et 1

  double gammaL = BO->eos[liq]->getGamma();
  double pInfL = BO->eos[liq]->getPInf();
  double cvL = BO->eos[liq]->getCv();
  double e0L = BO->eos[liq]->getERef();
  double s0L = BO->eos[liq]->getSRef();

  double gammaV = BO->eos[vap]->getGamma();
  double pInfV = BO->eos[vap]->getPInf();
  double cvV = BO->eos[vap]->getCv();
  double e0V = BO->eos[vap]->getERef();
  double s0V = BO->eos[vap]->getSRef();

  double A, B, C, D;
  A = (gammaL*cvL - gammaV*cvV + s0V - s0L) / (gammaV*cvV - cvV);
  B = (e0L - e0V) / (gammaV*cvV - cvV);
  C = (gammaV*cvV-gammaL*cvL)/ (gammaV*cvV - cvV);
  D = (gammaL*cvL - cvL) / (gammaV*cvV - cvV);

  //Processus iteratif recherche de la temperature de saturation
  int iteration(0);
  double Tsat(0.1*B/C);
  double f(0.), df(1.);
  do {
    Tsat -= f / df; iteration++;
    if (iteration > 50) {
      erreurs.push_back(Erreurs("nombre iterations trop grand dans recherche Tsat", __FILE__, __LINE__));
      break;
    }
    f = A+B/Tsat+C*log(Tsat)-log(pression +pInfV)+D*log(pression +pInfL);
    df = C/Tsat-B/(Tsat*Tsat);
  } while (fabs(f)>1e-10);

  double dfdp = -1. / (pression + pInfV) + D / (pression + pInfL);
  *dTsat = -dfdp / df;
  return Tsat;
}

//***********************************************************************

void FluxKapila::correctionEnergie(Cellule *cell, const int &nombrePhases, Prim type) const
{
  Phase *phase;

  //Extraction des donnees utiles
  for (int k = 0; k < nombrePhases; k++){
    phase = cell->getPhase(k, type);
    BO->ak[k] = phase->getAlpha();
    BO->rhok[k] = phase->getDensite();
  }
  //Calcul de le pression de melange via EOS de melange
  double rhoe = cell->getMelange(type)->getDensite() * cell->getMelange(type)->getEnergie();
  double p(rhoe), denom(0.), gamPinfSurGamMoinsUn(0.), eRef(0.), unSurGamMoinsUn(0.);
  for (int k = 0; k < nombrePhases; k++) {
    BO->eos[k]->renvoiSpecialEosMelange(gamPinfSurGamMoinsUn, eRef, unSurGamMoinsUn);
    p -= BO->ak[k]*(gamPinfSurGamMoinsUn + BO->rhok[k] * eRef);
    denom += BO->ak[k] * unSurGamMoinsUn;
  }
  p /= denom;
  //Remplacement des pressions par la pression calcule via EOS de melange
  for (int k = 0; k < nombrePhases; k++){
    phase = cell->getPhase(k, type);
    phase->setPression(p);
    phase->verifiePhase();
  }
  cell->getMelange(type)->setPression(p);
}

//***********************************************************************

void FluxKapila::prepareSourceAxi(Cellule *cell, const int &nombrePhases, double &r, double &v)
{
  sourceConsKap->miseAZero(nombrePhases);
  for (int k = 0;k<nombrePhases;k++)
  {
    sourceConsKap->m_alpha[k] = -m_alpha[k] * v / r;
    sourceConsKap->m_masse[k] = -m_masse[k] * v / r;
    sourceConsKap->m_energ[k] = -m_energ[k] * v / r;
  }
  sourceConsKap->m_qdm = -v / r * m_qdm;
  sourceConsKap->m_energMelange = -m_energMelange* v /r;
}

//***********************************************************************

void FluxKapila::prepareSourceGravite(Cellule *cell, const int &nombrePhases, Coord &Fg)
{
  sourceConsKap->miseAZero(nombrePhases);
  double rho = cell->getMelange()->getDensite();
  Coord u(m_qdm/rho);
  sourceConsKap->m_qdm = Fg;
  sourceConsKap->m_energMelange = Fg.scalaire(u);
}

//***********************************************************************

void FluxKapila::integreTermeSource(Cellule *cell, const double &dt, const int &nombrePhases)
{

  ////Version avec verification
  //for (int k = 0;k<nombrePhases;k++)
  //{
  //  if(1e-8 < m_alpha[k] + dt*sourceConsKap.m_alpha[k] < 1.0-1e-8) m_alpha[k] += dt*sourceConsKap.m_alpha[k];
  //  if(m_masse[k] + dt*sourceConsKap.m_masse[k] > 1e-8)  m_masse[k] += dt*sourceConsKap.m_masse[k];
  //  m_energ[k] += dt*sourceConsKap.m_energ[k];
  //}
  //m_qdm += dt*sourceConsKap.m_qdm;
  //m_energMelange += dt*sourceConsKap.m_energMelange;


  //Version simple mais non robuste
  for(int k=0;k<nombrePhases;k++)
  {
    m_alpha[k] += dt*sourceConsKap->m_alpha[k];
    m_masse[k] += dt*sourceConsKap->m_masse[k];
    m_energ[k] += dt*sourceConsKap->m_energ[k];
  }
  m_qdm += dt*sourceConsKap->m_qdm;
  m_energMelange += dt*sourceConsKap->m_energMelange;

}

//***********************************************************************

double FluxKapila::getAlpha(const int &numPhase) const
{
  return m_alpha[numPhase];
}

//***********************************************************************

double FluxKapila::getMasse(const int &numPhase) const
{
  return m_masse[numPhase];
}

//***********************************************************************

double FluxKapila::getEnergie(const int &numPhase) const
{
  return m_energ[numPhase];
}

//***********************************************************************

Coord FluxKapila::getQdm() const
{
  return m_qdm;
}

//***********************************************************************

double FluxKapila::getEnergieMel() const
{
	return m_energMelange;
}

//***********************************************************************

void FluxKapila::setCons(const Flux *cons, const int &nombrePhases)
{
  for (int k = 0; k < nombrePhases; k++) {
    m_alpha[k] = cons->getAlpha(k);
    m_masse[k] = cons->getMasse(k);
    m_energ[k] = cons->getEnergie(k);
  }
  m_qdm = cons->getQdm();
  m_energMelange = cons->getEnergieMel();
}

//***********************************************************************