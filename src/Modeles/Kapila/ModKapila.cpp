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
#include "ModKapila.h"
#include "PhaseKapila.h"

using namespace std;

const std::string ModKapila::NOM = "KAPILA";

//***********************************************************************

ModKapila::ModKapila(int &nombreTransports, const int &nombrePhases) :
  Modele(NOM,nombreTransports)
{
  fluxTempKapila = new FluxKapila(nombrePhases);
  sourceConsKap = new FluxKapila(nombrePhases);
  //fluxTempKapila->alloue(nombrePhases);
  //sourceConsKap.alloue(nombrePhases);
  m_rhokStar = new double[nombrePhases];
  m_pkStar = new double[nombrePhases];
  m_ekStar = new double[nombrePhases];
  m_EkStar = new double[nombrePhases];
  m_vkStar = new double[nombrePhases];
  m_YkStar = new double[nombrePhases];
  m_Hk0    = new double[nombrePhases];
  m_Yk0    = new double[nombrePhases];
}

//***********************************************************************

ModKapila::~ModKapila()
{
  delete fluxTempKapila;
  delete sourceConsKap;
  //fluxTempKapila->desalloue();
  //sourceConsKap.desalloue();
  delete[]m_rhokStar;
  delete[]m_pkStar;
  delete[]m_ekStar;
  delete[]m_EkStar;
  delete[]m_vkStar;
  delete[]m_YkStar;
  delete[]m_Hk0;
  delete[]m_Yk0;
}

//***********************************************************************

void ModKapila::alloueCons(Flux **cons, const int &nombrePhases)
{
  *cons = new FluxKapila(nombrePhases);
  //(*cons)->alloue(nombrePhases);
}

//***********************************************************************

void ModKapila::allouePhase(Phase **phase)
{
  *phase = new PhaseKapila;
}

//***********************************************************************

void ModKapila::alloueMelange(Melange **melange)
{
  *melange = new MelKapila;
}

//***********************************************************************

void ModKapila::alloueEos(Cellule &cell, const int &nombrePhases)
{
	for (int k = 0; k < nombrePhases; k++) { BO->eos[k] = cell.getPhase(k)->getEos(); }
}

//****************************************************************************
//************Probleme de Riemann entre deux mailles de calcul****************
//****************************************************************************

void ModKapila::resolRiemannInterne(Cellule &cellGauche, Cellule &cellDroite, const int &nombrePhases, const double &dxGauche, const double &dxDroite, double &dtMax) const
{
  Phase *vecPhase;
  double sL, sR;
  double pStar(0.), rhoStar(0.), uStar(0.), vStar(0.), wStar(0.), EStar(0.), eStar(0.);
  //Allocation des tableaux

  //Raccourcis
  double uL = cellGauche.getMelange()->getVitesse().getX(), cL = cellGauche.getMelange()->getVitesseSonFigee(), pL = cellGauche.getMelange()->getPression(), rhoL = cellGauche.getMelange()->getDensite();
  double uR = cellDroite.getMelange()->getVitesse().getX(), cR = cellDroite.getMelange()->getVitesseSonFigee(), pR = cellDroite.getMelange()->getPression(), rhoR = cellDroite.getMelange()->getDensite();

  //Davies
  sL = min(uL - cL, uR - cR);
  sR = max(uR + cR, uL + cL);

  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
  if (fabs(sR)>1.e-3) dtMax = min(dtMax, dxDroite / fabs(sR));

  //Calcul debits masse gauche et droite et sM
  double mL(rhoL*(sL - uL)), mR(rhoR*(sR - uR)), mkL, mkR;
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));
  if (fabs(sM)<1.e-8) sM = 0.;

  //Echantillonage solution
  if (sL >= 0.){
    for (int k = 0; k < nombrePhases; k++) {
      vecPhase = cellGauche.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double densite = vecPhase->getDensite();
      double energie = vecPhase->getEnergie();
      fluxTempKapila->m_alpha[k] = alpha*sM;
      fluxTempKapila->m_masse[k] = alpha*densite*uL;
      fluxTempKapila->m_energ[k] = alpha*densite*energie*uL;
    }
    double vitY = cellGauche.getMelange()->getVitesse().getY(); double vitZ = cellGauche.getMelange()->getVitesse().getZ();
    double energieTotale = cellGauche.getMelange()->getEnergieTotale();
    fluxTempKapila->m_qdm.setX(rhoL*uL*uL + pL);
    fluxTempKapila->m_qdm.setY(rhoL*vitY*uL);
    fluxTempKapila->m_qdm.setZ(rhoL*vitZ*uL);
    fluxTempKapila->m_energMelange = (rhoL*energieTotale + pL)*uL;

  }
  else if (sR <= 0.){
    for (int k = 0; k < nombrePhases; k++) {
      vecPhase = cellDroite.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double densite = vecPhase->getDensite();
      double energie = vecPhase->getEnergie();
      fluxTempKapila->m_alpha[k] = alpha*sM;
      fluxTempKapila->m_masse[k] = alpha*densite*uR;
      fluxTempKapila->m_energ[k] = alpha*densite*energie*uR;
    }
    double vitY = cellDroite.getMelange()->getVitesse().getY(); double vitZ = cellDroite.getMelange()->getVitesse().getZ();
    double energieTotale = cellDroite.getMelange()->getEnergieTotale();
    fluxTempKapila->m_qdm.setX(rhoR*uR*uR + pR);
    fluxTempKapila->m_qdm.setY(rhoR*vitY*uR);
    fluxTempKapila->m_qdm.setZ(rhoR*vitZ*uR);
    fluxTempKapila->m_energMelange = (rhoR*energieTotale + pR)*uR;

  }
  else if (sM >= 0.){
    //Calcul de l etat solution gauche
    double vitY = cellGauche.getMelange()->getVitesse().getY(); double vitZ = cellGauche.getMelange()->getVitesse().getZ();
    double energieTotale = cellGauche.getMelange()->getEnergieTotale();
    rhoStar = mL / (sL - sM);
    EStar = energieTotale + (sM - uL)*(sM + pL / mL);
    pStar = mL*(sM - uL) + pL;
    for (int k = 0; k < nombrePhases; k++) {
      vecPhase = cellGauche.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double densite = vecPhase->getDensite();
      double pression = vecPhase->getPression();
      energieTotale = vecPhase->getEnergieTotale();
      mkL = densite*(sL - uL);
      m_rhokStar[k] = mkL / (sL - sM);
      m_pkStar[k] = BO->eos[k]->calculPressionIsentrope(pression, densite, m_rhokStar[k]);
      m_ekStar[k] = BO->eos[k]->calculEnergie(m_rhokStar[k], m_pkStar[k]);
      m_EkStar[k] = energieTotale + (sM - uL)*(sM + pL / mkL); //non utilise (pour test sur Energies totales phases)
      fluxTempKapila->m_alpha[k] = alpha*sM;
      fluxTempKapila->m_masse[k] = alpha*m_rhokStar[k] * sM;
      fluxTempKapila->m_energ[k] = alpha*m_rhokStar[k] * m_ekStar[k] * sM;
    }
    fluxTempKapila->m_qdm.setX(rhoStar*sM*sM + pStar);
    fluxTempKapila->m_qdm.setY(rhoStar*vitY*sM);
    fluxTempKapila->m_qdm.setZ(rhoStar*vitZ*sM);
    fluxTempKapila->m_energMelange = (rhoStar*EStar + pStar)*sM;
  }
  else{
    //Calcul de l etat solution droite
    double vitY = cellDroite.getMelange()->getVitesse().getY(); double vitZ = cellDroite.getMelange()->getVitesse().getZ();
    double energieTotale = cellDroite.getMelange()->getEnergieTotale();
    rhoStar = mR / (sR - sM);
    EStar = energieTotale + (sM - uR)*(sM + pR / mR);
    pStar = mR*(sM - uR) + pR;
    for (int k = 0; k < nombrePhases; k++) {
      vecPhase = cellDroite.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double densite = vecPhase->getDensite();
      double pression = vecPhase->getPression();
      energieTotale = vecPhase->getEnergieTotale();
      mkR = densite*(sR - uR);
      m_rhokStar[k] = mkR / (sR - sM);
      m_pkStar[k] = BO->eos[k]->calculPressionIsentrope(pression, densite, m_rhokStar[k]);
      m_ekStar[k] = BO->eos[k]->calculEnergie(m_rhokStar[k], m_pkStar[k]);
      m_EkStar[k] = energieTotale + (sM - uR)*(sM + pR / mkR); //non utilise (pour test sur Energies totales phases)
      fluxTempKapila->m_alpha[k] = alpha*sM;
      fluxTempKapila->m_masse[k] = alpha*m_rhokStar[k] * sM;
      fluxTempKapila->m_energ[k] = alpha*m_rhokStar[k] * m_ekStar[k] * sM;
    }
    fluxTempKapila->m_qdm.setX(rhoStar*sM*sM + pStar);
    fluxTempKapila->m_qdm.setY(rhoStar*vitY*sM);
    fluxTempKapila->m_qdm.setZ(rhoStar*vitZ*sM);
    fluxTempKapila->m_energMelange = (rhoStar*EStar + pStar)*sM;
  }

  //vitesse bord
  fluxTempKapila->m_sM = sM;
}

//****************************************************************************
//**************Problemes de Riemann pour conditions aux limites**************
//****************************************************************************

void ModKapila::resolRiemannMur(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax) const
{
  double sL;
  double pStar(0.);

  //Raccourcis
  double uL = cellGauche.getMelange()->getVitesse().getX(), cL = cellGauche.getMelange()->getVitesseSonFigee(), pL = cellGauche.getMelange()->getPression(), rhoL = cellGauche.getMelange()->getDensite();

  sL = min(uL - cL, -uL - cL);
  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));

  pStar = rhoL*(uL - sL)*uL + pL;

  for (int k = 0; k < nombrePhases; k++)
  {
    fluxTempKapila->m_alpha[k] = 0.;
    fluxTempKapila->m_masse[k] = 0.;
    fluxTempKapila->m_energ[k] = 0.;
  }
  fluxTempKapila->m_qdm.setX(pStar);
  fluxTempKapila->m_qdm.setY(0.);
  fluxTempKapila->m_qdm.setZ(0.);
  fluxTempKapila->m_energMelange = 0.;

  //Non conservatif
  fluxTempKapila->m_sM = 0.;
}

//****************************************************************************

void ModKapila::resolRiemannInj(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double m0, const double *ak0, const double *rhok0, const double *pk0) const
{
  double sL, zL;
  double pStar(0.), uStar(0.), rhoStar(0.);

  //Raccourcis
  double uL = cellGauche.getMelange()->getVitesse().getX(), cL = cellGauche.getMelange()->getVitesseSonFigee(), pL = cellGauche.getMelange()->getPression(), rhoL = cellGauche.getMelange()->getDensite();
  double vL = cellGauche.getMelange()->getVitesse().getY(), wL = cellGauche.getMelange()->getVitesse().getZ();

  //Calcul enthalpie totale etat injection
  double rho0 = cellGauche.getMelange()->calculDensite(ak0, rhok0, nombrePhases);
  double u0 = m0 / rho0;
  for (int k = 0;k < nombrePhases;k++) {
    m_Hk0[k] = BO->eos[k]->calculEnthalpieTotale(rhok0[k], pk0[k], u0);
    m_Yk0[k] = ak0[k] * rhok0[k] / rho0;
  }

  //PORCEDE ITERATIF CALCUL DE LA PRESSION SOLUTION
  //-----------------------------------------------
  //Estimations acoustique onde sL
  sL = uL - cL;
  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
  zL = rhoL*cL;

  int iteration(0);
  pStar = pL;
  double f(0.), df(1.);
  double u, du, v, dv, hk;

  do {
    pStar -= f / df; iteration++;
    if (iteration > 50) Erreurs::messageErreur("nombre iterations trop grand dans resolRiemannInj modKapila");
    //Verification pression physique ?
    for (int k = 0; k < nombrePhases; k++) {
      BO->eos[k]->verifieEtCorrigePression(pStar);
    }
    //Relations acoustiques a gauche
    u = uL + (pL - pStar) / zL;
    if (u >= -1e-6) u = -1e-6;
    du = -1. / zL;
    //Calcul a partir de m0, Hk0, Yk0 a droite
    v = u / m0;
    dv = du / m0;
    f = v;
    df = dv;
    for (int k = 0; k < nombrePhases; k++) {
      hk = m_Hk0[k] - 0.5*u*u;
      m_vkStar[k] = BO->eos[k]->vfpfh(pStar, hk);
      double dvk = BO->eos[k]->dvdpch(pStar, hk) - BO->eos[k]->dvdhcp(pStar, hk)*u*du;
      f -= m_Yk0[k] * m_vkStar[k];
      df -= m_Yk0[k] * dvk;
    }
  } while (fabs(f)>1e-10);

  //On complete le Flux
  double Estar(0.5*(u*u + vL*vL + wL*wL)), ek, rhok;
  for (int k = 0; k<nombrePhases; k++) {
    rhok = 1. / m_vkStar[k];
    ek = BO->eos[k]->calculEnergie(rhok, pStar); Estar += m_Yk0[k] * ek;
    fluxTempKapila->m_alpha[k] = m_Yk0[k] * m_vkStar[k] / v*u;
    fluxTempKapila->m_masse[k] = fluxTempKapila->m_alpha[k] * rhok;
    fluxTempKapila->m_energ[k] = fluxTempKapila->m_alpha[k] * rhok*ek;
  }
  fluxTempKapila->m_qdm.setX(u*u / v + pStar);
  fluxTempKapila->m_qdm.setY(u*vL / v);
  fluxTempKapila->m_qdm.setZ(u*wL / v);
  fluxTempKapila->m_energMelange = (Estar / v + pStar)*u;

  //Non conservatif
  fluxTempKapila->m_sM = u;
}

//****************************************************************************

void ModKapila::resolRiemannRes(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double *ak0, const double *rhok0, const double *pk0) const
{
  double sL, zL;
  double pStar(0.), uStar(0.), rhoStar(0.), vStar(0.), wStar(0.);
  Phase *vecPhase;

  //Raccourcis
  double uL = cellGauche.getMelange()->getVitesse().getX(), cL = cellGauche.getMelange()->getVitesseSonWood(), pL = cellGauche.getMelange()->getPression(), rhoL = cellGauche.getMelange()->getDensite();
  double vL = cellGauche.getMelange()->getVitesse().getY(), wL = cellGauche.getMelange()->getVitesse().getZ();

  zL = rhoL*cL;

  //Estimation de la vitesse d onde sL en utilisant pStar = p0
  //----------------------------------------------------------
  double p0(cellGauche.getMelange()->calculPression(ak0, pk0, nombrePhases));
  pStar = p0;
  double v(0.), vmv0, mL, u;
  for (int k = 0; k < nombrePhases; k++) {
    vecPhase = cellGauche.getPhase(k);
    m_rhokStar[k] = BO->eos[k]->calculDensiteIsentrope(vecPhase->getPression(), vecPhase->getDensite(), pStar); //FP//DEV// mettre hugoniot
    //m_rhokStar[k] = BO->eos[k]->calculDensiteHugoniot(vecPhase->getPression(), vecPhase->getDensite(), pStar);
    v += vecPhase->getAlpha()*vecPhase->getDensite() / rhoL / m_rhokStar[k];
  }
  vmv0 = v - 1. / rhoL;
  if (fabs(vmv0) > 1e-10) {
    mL = sqrt((pL - p0) / vmv0);
  }
  else {
    mL = zL;
  }
  sL = uL - mL / rhoL;
  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
  u = uL + mL*vmv0;

  //Cas pathologiques
  //-----------------
  if (sL >= 0.) { //Cas sortie supersonique => etat gauche solution
    uStar = uL;
    pStar = pL;
    for (int k = 0; k < nombrePhases; k++) {
      vecPhase = cellGauche.getPhase(k);
      m_rhokStar[k] = vecPhase->getDensite();
      m_YkStar[k] = vecPhase->getAlpha()*vecPhase->getDensite() / rhoL;
    }
    rhoStar = rhoL;
    vStar = vL;
    wStar = wL;
  }
  else if (u >= -1e-3) { //Cas sortie subsonique => etat star gauche solution //FP//DEV//Ajouter Yk0 !
    uStar = u;
    pStar = p0;
    for (int k = 0; k < nombrePhases; k++) {
      vecPhase = cellGauche.getPhase(k);
      m_YkStar[k] = vecPhase->getAlpha()*vecPhase->getDensite() / rhoL;
    }
    rhoStar = 1. / v;
    vStar = vL;
    wStar = wL;
  }
  //Cas reservoir
  //-------------
  else { //Entree reservoir multiphasique => etat star droite solution

     //Calcul enthalpie totale etat reservoir
    double H0(0.), v0(0.);
    double rho0 = cellGauche.getMelange()->calculDensite(ak0, rhok0, nombrePhases);
    for (int k = 0;k < nombrePhases;k++) {
      m_Hk0[k] = BO->eos[k]->calculEnthalpieTotale(rhok0[k], pk0[k], v0); //vitesse nulle dans reservoir
      m_Yk0[k] = ak0[k] * rhok0[k] / rho0;
      H0 += m_Yk0[k] * m_Hk0[k];
    }

    //PORCEDE ITERATIF CALCUL DE LA PRESSION SOLUTION
    //-----------------------------------------------
    int iteration(0);
    double p(0.5*p0);
    double f(0.), df(1.);
    double hk, dhk, drhok;
    double uStarR(0.), duStarR(0.), uStarL(0.), duStarL(0.);
    do {
      p -= f / df; iteration++;
      if (iteration > 50) Erreurs::messageErreur("nombre iterations trop grand dans resolRiemannRes modKapila");
      //Verification pression physique ?
      for (int k = 0; k < nombrePhases; k++) { BO->eos[k]->verifieEtCorrigePression(p); }
      if (p > p0) { p = p0 - 1e-6; }

      //Relations reservoir a droite (Hk=cte et sk=cste)
      uStarR = H0; duStarR = 0.;
      for (int k = 0; k < nombrePhases; k++) {
        m_rhokStar[k] = BO->eos[k]->calculDensiteIsentrope(pk0[k], rhok0[k], p);
        hk = BO->eos[k]->calculEnthalpieIsentrope(pk0[k], rhok0[k], p, &dhk);
        uStarR -= m_Yk0[k] * hk;
        duStarR -= m_Yk0[k] * dhk;
      }
      uStarR = -sqrt(2.*uStarR);
      duStarR = duStarR / uStarR;

      //Relations a gauche (isentropiques ou Hugoniots)
      v = 0.; double dv(0.), dmL;
      double rhok;
      for (int k = 0; k < nombrePhases; k++) {
        vecPhase = cellGauche.getPhase(k);
        rhok = BO->eos[k]->calculDensiteIsentrope(vecPhase->getPression(), vecPhase->getDensite(), p, &drhok); //FP//DEV// mettre hugoniot
        //rhok = BO->eos[k]->calculDensiteHugoniot(vecPhase->getPression(), vecPhase->getDensite(), p, &drhok);
        double YkL = vecPhase->getAlpha()*vecPhase->getDensite() / rhoL;
        v += YkL / rhok;
        dv -= YkL / (rhok * rhok) * drhok;
      }
      vmv0 = v - 1. / rhoL;
      if (fabs(vmv0) > 1e-10) {
        mL = sqrt((pL - p) / vmv0);
        dmL = 0.5*(-vmv0 + (p - pL)*dv) / (vmv0*vmv0) / mL;
      }
      else {
        mL = zL;
        dmL = 0.;
      }
      sL = uL - mL / rhoL;
      if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
      uStarL = uL + mL*vmv0;
      duStarL = dmL*vmv0 + mL*dv;

      //Fonction a resoudre
      f = uStarR - uStarL;
      df = duStarR - duStarL;

    } while (fabs(f)>1e-3);

    pStar = p;
    uStar = 0.5*(uStarL + uStarR);
    rhoStar = 0.;
    for (int k = 0; k < nombrePhases; k++) { 
      m_YkStar[k] = m_Yk0[k];
      rhoStar += m_YkStar[k] / m_rhokStar[k]; 
    }
    rhoStar = 1. / rhoStar;
    vStar = 0.;
    wStar = 0.;
  }

  //On complete le Flux
  double EStar(0.5*(uStar*uStar + vStar*vStar + wStar*wStar)), ek; //FP//DEV// Vitesse transverse nulles dans le reservoir, non ?
  for (int k = 0; k < nombrePhases; k++) {
    ek = BO->eos[k]->calculEnergie(m_rhokStar[k], pStar); EStar += m_YkStar[k] * ek;
    fluxTempKapila->m_alpha[k] = m_YkStar[k] * rhoStar / m_rhokStar[k] * uStar;
    fluxTempKapila->m_masse[k] = rhoStar*m_YkStar[k] * uStar;
    fluxTempKapila->m_energ[k] = fluxTempKapila->m_masse[k] * ek;
  }
  fluxTempKapila->m_qdm.setX(rhoStar*uStar*uStar + pStar);
  fluxTempKapila->m_qdm.setY(rhoStar*uStar*vStar);
  fluxTempKapila->m_qdm.setZ(rhoStar*uStar*wStar);
  fluxTempKapila->m_energMelange = (rhoStar*EStar + pStar)*uStar;

  //Non conservatif
  fluxTempKapila->m_sM = uStar;
}

//****************************************************************************

void ModKapila::resolRiemannSortie(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double p0, double *debitSurf) const
{
  double sL, zL;
  double pStar(p0);
  Phase *vecPhase;

  //Raccourcis
  double uL = cellGauche.getMelange()->getVitesse().getX(), cL = cellGauche.getMelange()->getVitesseSonFigee(), pL = cellGauche.getMelange()->getPression(), rhoL = cellGauche.getMelange()->getDensite();
  double vL = cellGauche.getMelange()->getVitesse().getY(), wL = cellGauche.getMelange()->getVitesse().getZ();

  //Estimations acoustique onde sL
  //------------------------------
  zL = rhoL*cL;
  double v(0.), vmv0, mL, u;
  for (int k = 0; k < nombrePhases; k++) {
    vecPhase = cellGauche.getPhase(k);
    m_vkStar[k] = 1.0 / BO->eos[k]->calculDensiteIsentrope(vecPhase->getPression(), vecPhase->getDensite(), p0); //FP//DEV// mettre hugoniot
    v += vecPhase->getAlpha()*vecPhase->getDensite()/rhoL * m_vkStar[k];
  }
  vmv0 = v - 1. / rhoL;
  if (fabs(vmv0) > 1e-10) {
    mL = sqrt((pL - p0) / vmv0);
  }
  else {
    mL = zL;
  }
  sL = uL - mL / rhoL;
  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
  u = uL + mL*vmv0;

  //Cas patologique sL>0            //FP//Q// Peut etre cas u<0 a voir aussi
  if (sL >= 0.) { //Cas sortie supersonique => etat gauche solution
    u = uL;
    pStar = pL;
    for (int k = 0; k < nombrePhases; k++) m_vkStar[k] = 1. / cellGauche.getPhase(k)->getDensite();
    v = 1. / rhoL;
  }

  //On complete le Flux
  double Estar(0.5*(u*u + vL*vL + wL*wL)), ek, rhok;
  for (int k = 0; k < nombrePhases; k++) {
    vecPhase = cellGauche.getPhase(k);
    rhok = 1. / m_vkStar[k];
    double YkL = vecPhase->getAlpha()*vecPhase->getDensite() / rhoL;
    ek = BO->eos[k]->calculEnergie(rhok, pStar); Estar += YkL * ek;
    fluxTempKapila->m_alpha[k] = YkL * m_vkStar[k] / v * u;
    fluxTempKapila->m_masse[k] = fluxTempKapila->m_alpha[k] * rhok;
    fluxTempKapila->m_energ[k] = fluxTempKapila->m_alpha[k] * rhok*ek;
  }
  fluxTempKapila->m_qdm.setX(u*u / v + pStar);
  fluxTempKapila->m_qdm.setY(u*vL / v);
  fluxTempKapila->m_qdm.setZ(u*wL / v);
  fluxTempKapila->m_energMelange = (Estar / v + pStar)*u;

  //Non conservatif
  fluxTempKapila->m_sM = u;

  //Sauvegarde du debit surfacique
  for (int k = 0; k < nombrePhases; k++) {
    debitSurf[k] = fluxTempKapila->m_masse[k];
  }
}

//****************************************************************************
//***************Calcul des flux pour les equations de transport**************
//****************************************************************************

void ModKapila::resolRiemannTransportInterne(Cellule &cellGauche, Cellule &cellDroite, const int &nombreTransports)
{
	for (int k = 0; k < nombreTransports; k++) {
		fluxTempTransport[k].resolRiemann(cellGauche.getTransport(k).getValeur(), cellDroite.getTransport(k).getValeur(), fluxTempKapila->m_sM);
	}
}

//****************************************************************************

void ModKapila::resolRiemannTransportMur(const int &nombreTransports)
{
	for (int k = 0; k < nombreTransports; k++) {
		fluxTempTransport[k].resolRiemannMur();
	}
}

//****************************************************************************

void ModKapila::resolRiemannTransportInj(Cellule &cellGauche, const int &nombreTransports, double *valeurTransports)
{
	for (int k = 0; k < nombreTransports; k++) {
		fluxTempTransport[k].resolRiemannInj(cellGauche.getTransport(k).getValeur(), fluxTempKapila->m_sM, valeurTransports[k]);
	}
}

//****************************************************************************

void ModKapila::resolRiemannTransportRes(Cellule &cellGauche, const int &nombreTransports, double *valeurTransports)
{
	for (int k = 0; k < nombreTransports; k++) {
		fluxTempTransport[k].resolRiemannRes(cellGauche.getTransport(k).getValeur(), fluxTempKapila->m_sM, valeurTransports[k]);
	}
}

//****************************************************************************

void ModKapila::resolRiemannTransportSortie(Cellule &cellGauche, const int &nombreTransports, double *valeurTransports)
{
	for (int k = 0; k < nombreTransports; k++) {
		fluxTempTransport[k].resolRiemannSortie(cellGauche.getTransport(k).getValeur(), fluxTempKapila->m_sM, valeurTransports[k]);
	}
}

//****************************************************************************

double ModKapila::getSM()
{
  return fluxTempKapila->m_sM;
}

//****************************************************************************
//*****************************Autres methodes********************************
//****************************************************************************

void ModKapila::projectionRepereAbsolu(const Coord normale, const Coord tangente, const Coord binormale) const
{
  Coord fluxProjete;
  fluxProjete.setX(normale.getX()*fluxTempKapila->m_qdm.getX() + tangente.getX()*fluxTempKapila->m_qdm.getY() + binormale.getX()*fluxTempKapila->m_qdm.getZ());
  fluxProjete.setY(normale.getY()*fluxTempKapila->m_qdm.getX() + tangente.getY()*fluxTempKapila->m_qdm.getY() + binormale.getY()*fluxTempKapila->m_qdm.getZ());
  fluxProjete.setZ(normale.getZ()*fluxTempKapila->m_qdm.getX() + tangente.getZ()*fluxTempKapila->m_qdm.getY() + binormale.getZ()*fluxTempKapila->m_qdm.getZ());
  fluxTempKapila->m_qdm.setXYZ(fluxProjete.getX(), fluxProjete.getY(), fluxProjete.getZ());
}

//****************************************************************************

string ModKapila::quiSuisJe() const
{
  return m_nom;
}

//****************************************************************************