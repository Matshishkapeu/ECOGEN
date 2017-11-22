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
#include <string>
#include "ModEuler.h"
#include "PhaseEuler.h"

using namespace std;

const std::string ModEuler::NOM = "EULER";

//****************************************************************************

ModEuler::ModEuler(const int &nombreTransports) :
  Modele(NOM, nombreTransports)
{}

//****************************************************************************

ModEuler::~ModEuler(){}

//****************************************************************************

void ModEuler::alloueCons(Flux **cons, const int &nombrePhases)
{
  *cons = new FluxEuler;
}

//***********************************************************************

void ModEuler::allouePhase(Phase **phase)
{
  *phase = new PhaseEuler;
}

//***********************************************************************

void ModEuler::alloueMelange(Melange **melange)
{
  *melange = new MelEuler;
}

//****************************************************************************
//************Probleme de Riemann entre deux mailles de calcul****************
//****************************************************************************

void ModEuler::resolRiemannInterne(Cellule &cellGauche, Cellule &cellDroite, const int &nombrePhases, const double &dxGauche, const double &dxDroite, double &dtMax) const
{
  Eos *eos;

  double cL, cR, sL, sR;
  double uL, uR, vL, vR, wL, wR, pL, pR, rhoL, rhoR, EL, ER;

  Phase *phaseGauche(0), *phaseDroite(0);
  phaseGauche = cellGauche.getPhase(0);
  phaseDroite = cellDroite.getPhase(0);

  eos = phaseGauche->getEos();
  uL = phaseGauche->getU(); vL = phaseGauche->getV(); wL = phaseGauche->getW();
  pL = phaseGauche->getPression();
  rhoL = phaseGauche->getDensite();
  cL = phaseGauche->getVitesseSon();
  EL = phaseGauche->getEnergieTotale();

  eos = phaseDroite->getEos();
  uR = phaseDroite->getU(); vR = phaseDroite->getV(); wR = phaseDroite->getW();
  pR = phaseDroite->getPression();
  rhoR = phaseDroite->getDensite();
  cR = phaseDroite->getVitesseSon();
  ER = phaseDroite->getEnergieTotale();

  sL = min(uL - cL, uR - cR);
  sR = max(uR + cR, uL + cL);

  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
  if (fabs(sR)>1.e-3) dtMax = min(dtMax, dxDroite / fabs(sR));

  //Calcul debits masse gauche et droite et sM
  double mL(rhoL*(sL - uL)), mR(rhoR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));
  if (fabs(sM)<1.e-8) sM = 0.;

  if (sL > 0.){
    fluxTempEuler.m_masse = rhoL*uL;
    fluxTempEuler.m_qdm.setX(rhoL*uL*uL + pL);
    fluxTempEuler.m_qdm.setY(rhoL*vL*uL);
    fluxTempEuler.m_qdm.setZ(rhoL*wL*uL);
    fluxTempEuler.m_energ = (rhoL*EL + pL)*uL;
  }
  else if (sR < 0.){
    fluxTempEuler.m_masse = rhoR*uR;
    fluxTempEuler.m_qdm.setX(rhoR*uR*uR + pR);
    fluxTempEuler.m_qdm.setY(rhoR*vR*uR);
    fluxTempEuler.m_qdm.setZ(rhoR*wR*uR);
    fluxTempEuler.m_energ = (rhoR*ER + pR)*uR;
  }

  ////Option HLL
  //else if (fabs(sR - sL)>1.e-3)
  //{
  //  fluxTempEuler.m_masse = (rhoR*uR*sL - rhoL*uL*sR + sL*sR*(rhoL - rhoR)) / (sL - sR);
  //  fluxTempEuler.m_qdm.setX(((rhoR*uR*uR + pR)*sL - (rhoL*uL*uL + pL)*sR + sL*sR*(rhoL*uL - rhoR*uR)) / (sL - sR));
  //  fluxTempEuler.m_qdm.setY((rhoR*uR*vR*sL - rhoL*uL*vL*sR + sL*sR*(rhoL*vL - rhoR*vR)) / (sL - sR));
  //  fluxTempEuler.m_qdm.setZ((rhoR*uR*wR*sL - rhoL*uL*wL*sR + sL*sR*(rhoL*wL - rhoR*wR)) / (sL - sR));
  //  fluxTempEuler.m_energ = ((rhoR*ER + pR)*uR*sL - (rhoL*EL + pL)*uL*sR + sL*sR*(rhoL*EL - rhoR*ER)) / (sL - sR);
  //}

  //Option HLLC
  else if (sM >= 0.) {
    double pStar = mL*(sM - uL) + pL;
    double rhoStar = mL / (sL - sM);
    double Estar = EL + (sM - uL)*(sM + pL / mL);
    fluxTempEuler.m_masse = rhoStar*sM;
    fluxTempEuler.m_qdm.setX(rhoStar*sM*sM+pStar);
    fluxTempEuler.m_qdm.setY(rhoStar*sM*vL);
    fluxTempEuler.m_qdm.setZ(rhoStar*sM*wL);
    fluxTempEuler.m_energ = (rhoStar*Estar + pStar)*sM;
  }
  else {
    double pStar = mR*(sM - uR) + pR;
    double rhoStar = mR / (sR - sM);
    double Estar = ER + (sM - uR)*(sM + pR / mR);
    fluxTempEuler.m_masse = rhoStar*sM;
    fluxTempEuler.m_qdm.setX(rhoStar*sM*sM + pStar);
    fluxTempEuler.m_qdm.setY(rhoStar*sM*vR);
    fluxTempEuler.m_qdm.setZ(rhoStar*sM*wR);
    fluxTempEuler.m_energ = (rhoStar*Estar + pStar)*sM;
  }

  //vitesse bord
  fluxTempEuler.m_sM = sM;
}

//****************************************************************************
//**************Problemes de Riemann pour conditions aux limites**************
//****************************************************************************

void ModEuler::resolRiemannMur(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax) const
{
  Eos *eos;

  double cL, sL;
  double uL, pL, rhoL;
  double pStar(0.);

  Phase *phaseGauche(0);
  phaseGauche = cellGauche.getPhase(0);

  eos = phaseGauche->getEos();
  uL = phaseGauche->getU();
  pL = phaseGauche->getPression();
  rhoL = phaseGauche->getDensite();
  cL = phaseGauche->getVitesseSon();

  sL = min(uL - cL, -uL - cL);
  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));

  pStar = rhoL*uL*(uL - sL) + pL;

  fluxTempEuler.m_masse = 0.;
  fluxTempEuler.m_qdm.setX(pStar);
  fluxTempEuler.m_qdm.setY(0.);
  fluxTempEuler.m_qdm.setZ(0.);
  fluxTempEuler.m_energ = 0.;

  //vitesse bord
  fluxTempEuler.m_sM = 0.;
}

//****************************************************************************

void ModEuler::resolRiemannInj(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double m0, const double *ak0, const double *rhok0, const double *pk0) const
{
  Eos *eos;

  double H0, u0;

  double cL, sL, zL;
  double uL, pL, rhoL, vL, wL;
  double uStar(0.), rhoStar(0.), pStar(0.), eStar(0.);

  Phase *phaseGauche(0);
  phaseGauche = cellGauche.getPhase(0);

  eos = phaseGauche->getEos();
  uL = phaseGauche->getU();
  vL = phaseGauche->getV();
  wL = phaseGauche->getW();
  pL = phaseGauche->getPression();
  rhoL = phaseGauche->getDensite();
  cL = phaseGauche->getVitesseSon();
  zL = rhoL*cL;

  sL = uL - cL;
  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));

  u0 = m0 / rhok0[0];
  H0 = eos->calculEnthalpieTotale(rhok0[0], pk0[0], u0);
  
  //Version EOS generale
  //--------------------
  int iteration(0);
  double p(pL);
  double f(0.), df(1.);
  double u, du, v, dv, h;

  do{
    p -= f / df; iteration++;
    if (iteration > 50) Erreurs::messageErreur("nombre iterations trop grand dans resolRiemannInj modEuler");
    //Verification pression physique ?
    eos->verifieEtCorrigePression(p);
    //Relations acoustique a gauche
    u = uL + (pL - p) / zL;
    if (u >= -1e-6) u = -1e-6;
    du = -1. / zL;
    v = u / m0;
    dv = du / m0;
    f = v;
    df = dv;
    //Calcul a partir de m0, H0 a droite
    h = H0 - 0.5*u*u;
    v = eos->vfpfh(p, h);
    dv = eos->dvdpch(p, h) - eos->dvdhcp(p, h)*u*du;
    f -= v;
    df -= dv;

  } while (fabs(f)>1e-10);
  uStar = u;
  pStar = p;
  rhoStar = m0 / uStar;
  eStar = eos->calculEnergie(rhoStar, pStar);

  ////Version GP ou SG uniquement => calcul exact
  ////-------------------------------------------
  //double *dataEos;
  //double a, b, c, delta, u1, u2, gammaTemp;
  //eos->renvoiInfo(dataEos);
  //gammaTemp = (dataEos[0] - 1)* m0 / dataEos[0];
  //a = 0.5*gammaTemp - zL;
  //b = pL + zL*uL;
  //c = -gammaTemp*H0;
  //delta = b*b - 4 * a*c;
  //u1 = (-b - sqrt(delta)) / (2 * a);
  //u2 = (-b + sqrt(delta)) / (2 * a);
  //uStar = min(u1, u2);
  //pStar = pL + zL*uL - zL*uStar;
  //rhoStar = m0 / uStar;
  //eStar = eos->calculEnergie(rhoStar, pStar);

  fluxTempEuler.m_masse = rhoStar*uStar;
  fluxTempEuler.m_qdm.setX(rhoStar*uStar*uStar + pStar);
  fluxTempEuler.m_qdm.setY(rhoStar*uStar*vL);
  fluxTempEuler.m_qdm.setZ(rhoStar*uStar*wL);
  fluxTempEuler.m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vL*vL + wL*wL)) + pStar)*uStar;

  //vitesse bord
  fluxTempEuler.m_sM = uStar;
}

//****************************************************************************

void ModEuler::resolRiemannRes(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double *ak0, const double *rhok0, const double *pk0) const
{
  Eos *eos;

  double cL, sL, zL;
  double uL, pL, rhoL, vL, wL;
  double uStar(0.), rhoStar(0.), pStar(0.), eStar(0.), vStar(0.), wStar(0.);

  Phase *phaseGauche(0);
  phaseGauche = cellGauche.getPhase(0);

  eos = phaseGauche->getEos();
  uL = phaseGauche->getU();
  vL = phaseGauche->getV();
  wL = phaseGauche->getW();
  pL = phaseGauche->getPression();
  rhoL = phaseGauche->getDensite();
  cL = phaseGauche->getVitesseSon();
  zL = rhoL*cL;

  //Estimation vitesse onde a gauche avec pStar = p0
  //------------------------------------------------
  pStar = pk0[0];
  double v(0.), vmv0, mL, u;
  v = 1./eos->calculDensiteIsentrope(pL, rhoL, pStar); //FP//DEV// mettre hugoniot
  vmv0 = v - 1. / rhoL;
  if (fabs(vmv0) > 1e-10) {
    mL = sqrt((pL - pk0[0]) / vmv0);
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
    rhoStar = rhoL;
    vStar = vL;
    wStar = wL;
  }
  else if (u >= -1e-3) { //Cas sortie subsonique => etat star gauche solution
    uStar = u;
    pStar = pk0[0];
    rhoStar = 1. / v;
    vStar = vL;
    wStar = wL;
  }
  //Cas reservoir
  //-------------
  else { //Entree reservoir => etat star droite solution
    //Calcul enthalpie totale etat reservoir
    double H0(0.);
    double v0 = 0.;
    H0 = eos->calculEnthalpieTotale(rhok0[0], pk0[0], v0);

    //PORCEDE ITERATIF CALCUL DE LA PRESSION SOLUTION
    //-----------------------------------------------
    int iteration(0);
    double p(0.5*pk0[0]);
    double f(0.), df(1.);
    double dv, h, dh, drho;
    double uStarR(0.), duStarR(0.), uStarL(0.), duStarL(0.);
    do {
      p -= f / df; iteration++;
      if (iteration > 50) Erreurs::messageErreur("nombre iterations trop grand dans resolRiemannRes modEuler");
      //Verification pression physique ?
      eos->verifieEtCorrigePression(p);
      if (p > pk0[0]) { p = pk0[0] - 1e-6; }

      //Relations reservoir a droite (H=cte et s=cste)
      rhoStar = eos->calculDensiteIsentrope(pk0[0], rhok0[0], p);
      h = eos->calculEnthalpieIsentrope(pk0[0], rhok0[0], p, &dh);
      uStarR = -sqrt(2.*(H0 - h));
      duStarR = -dh / uStarR; ;

      //Relations a gauche (isentropiques //FP//DEV// chocs a mettre)
      double dmL;
      v = 1.0 / eos->calculDensiteIsentrope(pL, rhoL, p, &drho); //FP//DEV// mettre hugoniot
      //v = 1.0 / eos->calculDensiteHugoniot(pL, rhoL, p, &drho);
      dv = - v*v*drho;
      vmv0 = v - 1. / rhoL;
      if (fabs(vmv0) > 1e-10) {
        mL = sqrt((pL - p) / vmv0);
        dmL = 0.5*(-vmv0 + (p - pL)*dv) / (vmv0*vmv0) / mL;
      }
      else { //si depassement densite limite sous choc => relations acoustiques
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
    vStar = 0.;
    wStar = 0.;
  }

  eStar = eos->calculEnergie(rhoStar, pStar);

  fluxTempEuler.m_masse = rhoStar*uStar;
  fluxTempEuler.m_qdm.setX(rhoStar*uStar*uStar + pStar);
  fluxTempEuler.m_qdm.setY(rhoStar*uStar*vStar);
  fluxTempEuler.m_qdm.setZ(rhoStar*uStar*wStar);
  fluxTempEuler.m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vStar*vStar + wStar*wStar)) + pStar)*uStar; //FP//DEV// Vitesse transverse nulles dans le reservoir, non ?

  //vitesse bord
  fluxTempEuler.m_sM = uStar;
}

//****************************************************************************

void ModEuler::resolRiemannSortie(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double p0, double *debitSurf) const
{
  Eos *eos;

  double cL, sL, zL;
  double uL, pL, rhoL, vL, wL;
  double uStar(0.), rhoStar(0.), pStar(0.), eStar(0.);

  double *dataEos;

  Phase *phaseGauche(0);
  phaseGauche = cellGauche.getPhase(0);

  eos = phaseGauche->getEos();
  uL = phaseGauche->getU();
  vL = phaseGauche->getV();
  wL = phaseGauche->getW();
  pL = phaseGauche->getPression();
  rhoL = phaseGauche->getDensite();
  cL = phaseGauche->getVitesseSon();
  zL = rhoL*cL;

  sL = uL - cL;
  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));

  eos->renvoiInfo(dataEos); //FP//TODO// A enlever : tres tres moche + probleme de desallocation !
  pStar = p0;
  rhoStar = rhoL*pow((pStar / pL), 1. / dataEos[0]);
  uStar = (pL + zL*uL - pStar) / zL;

  //Cas patologique sL>0            //FP//Q// Peut etre cas u<0 a voir aussi
  if (sL >= 0.) { //Cas sortie supersonique => etat gauche solution
    uStar = uL;
    pStar = pL;
    rhoStar = rhoL;
  }

  eStar = eos->calculEnergie(rhoStar, pStar);
  fluxTempEuler.m_masse = rhoStar*uStar;
  fluxTempEuler.m_qdm.setX(rhoStar*uStar*uStar + pStar);
  fluxTempEuler.m_qdm.setY(rhoStar*uStar*vL);
  fluxTempEuler.m_qdm.setZ(rhoStar*uStar*wL);
  fluxTempEuler.m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vL*vL + wL*wL)) + pStar)*uStar;

  //vitesse bord
  fluxTempEuler.m_sM = uStar;

  //Sauvegarde du debit surfacique
  debitSurf[0] = fluxTempEuler.m_masse;
}

//****************************************************************************

double ModEuler::getSM()
{
  return fluxTempEuler.m_sM;
}

//****************************************************************************

void ModEuler::projectionRepereAbsolu(const Coord normale, const Coord tangente, const Coord binormale) const
{
  Coord fluxProjete;
  fluxProjete.setX(normale.getX()*fluxTempEuler.m_qdm.getX() + tangente.getX()*fluxTempEuler.m_qdm.getY() + binormale.getX()*fluxTempEuler.m_qdm.getZ());
  fluxProjete.setY(normale.getY()*fluxTempEuler.m_qdm.getX() + tangente.getY()*fluxTempEuler.m_qdm.getY() + binormale.getY()*fluxTempEuler.m_qdm.getZ());
  fluxProjete.setZ(normale.getZ()*fluxTempEuler.m_qdm.getX() + tangente.getZ()*fluxTempEuler.m_qdm.getY() + binormale.getZ()*fluxTempEuler.m_qdm.getZ());
  fluxTempEuler.m_qdm.setXYZ(fluxProjete.getX(), fluxProjete.getY(), fluxProjete.getZ());
}

//****************************************************************************

string ModEuler::quiSuisJe() const
{
  return m_nom;
}

//****************************************************************************