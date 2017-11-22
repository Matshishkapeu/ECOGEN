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
#include "ModEulerHomogene.h"
#include "PhaseEulerHomogene.h"

using namespace std;

const std::string ModEulerHomogene::NOM = "EULERHOMOGENE";

//****************************************************************************

ModEulerHomogene::ModEulerHomogene(const int &nombreTransports, const int liquide, const int vapeur) :
  Modele(NOM, nombreTransports), m_liq(liquide-1), m_vap(vapeur-1)
{}

//****************************************************************************

ModEulerHomogene::~ModEulerHomogene(){}

//****************************************************************************

void ModEulerHomogene::alloueCons(Flux **cons, const int &nombrePhases)
{
  *cons = new FluxEulerHomogene(this);
}

//***********************************************************************

void ModEulerHomogene::allouePhase(Phase **phase)
{
  *phase = new PhaseEulerHomogene;
}

//***********************************************************************

void ModEulerHomogene::alloueMelange(Melange **melange)
{
  *melange = new MelEulerHomogene;
}

//****************************************************************************
//************Probleme de Riemann entre deux mailles de calcul****************
//****************************************************************************

void ModEulerHomogene::resolRiemannInterne(Cellule &cellGauche, Cellule &cellDroite, const int &nombrePhases, const double &dxGauche, const double &dxDroite, double &dtMax) const
{
  //Eos *eos;
  //Phase *vecPhase;
  double sL, sR;
  
  //Raccourcis //FP//TODO// Changer vitesse du son figee en vitesse du son equilibre !
  double uL = cellGauche.getMelange()->getVitesse().getX(), vL = cellGauche.getMelange()->getVitesse().getY(), wL = cellGauche.getMelange()->getVitesse().getZ(), cL = cellGauche.getMelange()->getVitesseSonFigee(), pL = cellGauche.getMelange()->getPression(), rhoL = cellGauche.getMelange()->getDensite();
  double uR = cellDroite.getMelange()->getVitesse().getX(), vR = cellDroite.getMelange()->getVitesse().getY(), wR = cellDroite.getMelange()->getVitesse().getZ(), cR = cellDroite.getMelange()->getVitesseSonFigee(), pR = cellDroite.getMelange()->getPression(), rhoR = cellDroite.getMelange()->getDensite();
  double EL = cellGauche.getMelange()->getEnergieTotale(), ER = cellDroite.getMelange()->getEnergieTotale();

  //Davies
  sL = min(uL - cL, uR - cR);
  sR = max(uR + cR, uL + cL);
  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
  if (fabs(sR)>1.e-3) dtMax = min(dtMax, dxDroite / fabs(sR));

  //sL = -1000.; //FP//TEST// 
  //sR = 1000.;

  ////Calcul des grandeurs melange gauche et droite
  //for (int k = 0; k < 2; k++) {
  //  rhoL += cellGauche.getPhase(k)->getAlpha()*cellGauche.getPhase(k)->getDensite();
  //  rhoR += cellDroite.getPhase(k)->getAlpha()*cellDroite.getPhase(k)->getDensite();

  //}
  
  //Calcul debits masse gauche et droite et sM
  double mL(rhoL*(sL - uL)), mR(rhoR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));
  if (fabs(sM)<1.e-8) sM = 0.;

  if (sL > 0.){
    fluxTempEulerHomogene.m_masse = rhoL*uL;
    fluxTempEulerHomogene.m_qdm.setX(rhoL*uL*uL + pL);
    fluxTempEulerHomogene.m_qdm.setY(rhoL*vL*uL);
    fluxTempEulerHomogene.m_qdm.setZ(rhoL*wL*uL);
    fluxTempEulerHomogene.m_energ = (rhoL*EL + pL)*uL;
  }
  else if (sR < 0.){
    fluxTempEulerHomogene.m_masse = rhoR*uR;
    fluxTempEulerHomogene.m_qdm.setX(rhoR*uR*uR + pR);
    fluxTempEulerHomogene.m_qdm.setY(rhoR*vR*uR);
    fluxTempEulerHomogene.m_qdm.setZ(rhoR*wR*uR);
    fluxTempEulerHomogene.m_energ = (rhoR*ER + pR)*uR;
  }

  ////Option HLL
  //else if (fabs(sR - sL)>1.e-3)
  //{
  //  fluxTempEulerHomogene.m_masse = (rhoR*uR*sL - rhoL*uL*sR + sL*sR*(rhoL - rhoR)) / (sL - sR);
  //  fluxTempEulerHomogene.m_qdm.setX(((rhoR*uR*uR + pR)*sL - (rhoL*uL*uL + pL)*sR + sL*sR*(rhoL*uL - rhoR*uR)) / (sL - sR));
  //  fluxTempEulerHomogene.m_qdm.setY((rhoR*uR*vR*sL - rhoL*uL*vL*sR + sL*sR*(rhoL*vL - rhoR*vR)) / (sL - sR));
  //  fluxTempEulerHomogene.m_qdm.setZ((rhoR*uR*wR*sL - rhoL*uL*wL*sR + sL*sR*(rhoL*wL - rhoR*wR)) / (sL - sR));
  //  fluxTempEulerHomogene.m_energ = ((rhoR*ER + pR)*uR*sL - (rhoL*EL + pL)*uL*sR + sL*sR*(rhoL*EL - rhoR*ER)) / (sL - sR);
  //}

  //Option HLLC
  else if (sM >= 0.) {
    double pStar = mL*(sM - uL) + pL;
    double rhoStar = mL / (sL - sM);
    double Estar = EL + (sM - uL)*(sM + pL / mL);
    fluxTempEulerHomogene.m_masse = rhoStar*sM;
    fluxTempEulerHomogene.m_qdm.setX(rhoStar*sM*sM+pStar);
    fluxTempEulerHomogene.m_qdm.setY(rhoStar*sM*vL);
    fluxTempEulerHomogene.m_qdm.setZ(rhoStar*sM*wL);
    fluxTempEulerHomogene.m_energ = (rhoStar*Estar + pStar)*sM;
  }
  else {
    double pStar = mR*(sM - uR) + pR;
    double rhoStar = mR / (sR - sM);
    double Estar = ER + (sM - uR)*(sM + pR / mR);
    fluxTempEulerHomogene.m_masse = rhoStar*sM;
    fluxTempEulerHomogene.m_qdm.setX(rhoStar*sM*sM + pStar);
    fluxTempEulerHomogene.m_qdm.setY(rhoStar*sM*vR);
    fluxTempEulerHomogene.m_qdm.setZ(rhoStar*sM*wR);
    fluxTempEulerHomogene.m_energ = (rhoStar*Estar + pStar)*sM;
  }

  //vitesse bord
  fluxTempEulerHomogene.m_sM = sM;
}

//****************************************************************************
//**************Problemes de Riemann pour conditions aux limites**************
//****************************************************************************

//void ModEulerHomogene::resolRiemannMur(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax) const
//{
//  Eos *eos;
//
//  double cL, sL;
//  double uL, pL, rhoL;
//  double pStar(0.);
//
//  Phase *phaseGauche(0);
//  phaseGauche = cellGauche.getPhase(0);
//
//  eos = phaseGauche->getEos();
//  uL = phaseGauche->getU();
//  pL = phaseGauche->getPression();
//  rhoL = phaseGauche->getDensite();
//  cL = eos->calculVitesseSon(rhoL, pL);
//
//  sL = min(uL - cL, -uL - cL);
//  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
//
//  pStar = rhoL*uL*(uL - sL) + pL;
//
//  fluxTempEulerHomogene.m_masse = 0.;
//  fluxTempEulerHomogene.m_qdm.setX(pStar);
//  fluxTempEulerHomogene.m_qdm.setY(0.);
//  fluxTempEulerHomogene.m_qdm.setZ(0.);
//  fluxTempEulerHomogene.m_energ = 0.;
//
//  //vitesse bord
//  fluxTempEulerHomogene.m_sM = 0.;
//}
//
////****************************************************************************
//
//void ModEulerHomogene::resolRiemannInj(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double m0, const double *ak0, const double *rhok0, const double *pk0) const
//{
//  Eos *eos;
//
//  double H0, u0;
//
//  double cL, sL, zL;
//  double uL, pL, rhoL, vL, wL;
//  double uStar(0.), rhoStar(0.), pStar(0.), eStar(0.);
//
//  Phase *phaseGauche(0);
//  phaseGauche = cellGauche.getPhase(0);
//
//  eos = phaseGauche->getEos();
//  uL = phaseGauche->getU();
//  vL = phaseGauche->getV();
//  wL = phaseGauche->getW();
//  pL = phaseGauche->getPression();
//  rhoL = phaseGauche->getDensite();
//  cL = eos->calculVitesseSon(rhoL, pL);
//  zL = rhoL*cL;
//
//  sL = uL - cL;
//  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
//
//  u0 = m0 / rhok0[0];
//  H0 = eos->calculEnthalpieTotale(rhok0[0], pk0[0], u0);
//  
//  //Version EOS generale
//  //--------------------
//  int iteration(0);
//  double p(pL);
//  double f(0.), df(1.);
//  double u, du, v, dv, h;
//
//  do{
//    p -= f / df; iteration++;
//    if (iteration > 50) Erreurs::messageErreur("nombre iterations trop grand dans resolRiemannInj modEuler");
//    //Verification pression physique ?
//    eos->verifieEtCorrigePression(p);
//    //Relations acoustique a gauche
//    u = uL + (pL - p) / zL;
//    if (u >= -1e-6) u = -1e-6;
//    du = -1. / zL;
//    v = u / m0;
//    dv = du / m0;
//    f = v;
//    df = dv;
//    //Calcul a partir de m0, H0 a droite
//    h = H0 - 0.5*u*u;
//    v = eos->vfpfh(p, h);
//    dv = eos->dvdpch(p, h) - eos->dvdhcp(p, h)*u*du;
//    f -= v;
//    df -= dv;
//
//  } while (fabs(f)>1e-10);
//  uStar = u;
//  pStar = p;
//  rhoStar = m0 / uStar;
//  eStar = eos->calculEnergie(rhoStar, pStar);
//
//  ////Version GP ou SG uniquement => calcul exact
//  ////-------------------------------------------
//  //double *dataEos;
//  //double a, b, c, delta, u1, u2, gammaTemp;
//  //eos->renvoiInfo(dataEos);
//  //gammaTemp = (dataEos[0] - 1)* m0 / dataEos[0];
//  //a = 0.5*gammaTemp - zL;
//  //b = pL + zL*uL;
//  //c = -gammaTemp*H0;
//  //delta = b*b - 4 * a*c;
//  //u1 = (-b - sqrt(delta)) / (2 * a);
//  //u2 = (-b + sqrt(delta)) / (2 * a);
//  //uStar = min(u1, u2);
//  //pStar = pL + zL*uL - zL*uStar;
//  //rhoStar = m0 / uStar;
//  //eStar = eos->calculEnergie(rhoStar, pStar);
//
//  fluxTempEulerHomogene.m_masse = rhoStar*uStar;
//  fluxTempEulerHomogene.m_qdm.setX(rhoStar*uStar*uStar + pStar);
//  fluxTempEulerHomogene.m_qdm.setY(rhoStar*uStar*vL);
//  fluxTempEulerHomogene.m_qdm.setZ(rhoStar*uStar*wL);
//  fluxTempEulerHomogene.m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vL*vL + wL*wL)) + pStar)*uStar;
//
//  //vitesse bord
//  fluxTempEulerHomogene.m_sM = uStar;
//}
//
////****************************************************************************
//
//void ModEulerHomogene::resolRiemannRes(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double *ak0, const double *rhok0, const double *pk0) const
//{
//  Eos *eos;
//
//  double cL, sL, zL;
//  double uL, pL, rhoL, vL, wL;
//  double uStar(0.), rhoStar(0.), pStar(0.), eStar(0.), vStar(0.), wStar(0.);
//
//  Phase *phaseGauche(0);
//  phaseGauche = cellGauche.getPhase(0);
//
//  eos = phaseGauche->getEos();
//  uL = phaseGauche->getU();
//  vL = phaseGauche->getV();
//  wL = phaseGauche->getW();
//  pL = phaseGauche->getPression();
//  rhoL = phaseGauche->getDensite();
//  cL = eos->calculVitesseSon(rhoL, pL);
//  zL = rhoL*cL;
//
//  //Estimation vitesse onde a gauche avec pStar = p0
//  //------------------------------------------------
//  pStar = pk0[0];
//  double v(0.), vmv0, mL, u;
//  v = 1./eos->calculDensiteIsentrope(pL, rhoL, pStar); //FP//DEV// mettre hugoniot
//  vmv0 = v - 1. / rhoL;
//  if (fabs(vmv0) > 1e-10) {
//    mL = sqrt((pL - pk0[0]) / vmv0);
//  }
//  else {
//    mL = zL;
//  }
//  sL = uL - mL / rhoL;
//  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
//  u = uL + mL*vmv0;
//
//  //Cas pathologiques
//  //-----------------
//  if (sL >= 0.) { //Cas sortie supersonique => etat gauche solution
//    uStar = uL;
//    pStar = pL;
//    rhoStar = rhoL;
//    vStar = vL;
//    wStar = wL;
//  }
//  else if (u >= -1e-3) { //Cas sortie subsonique => etat star gauche solution
//    uStar = u;
//    pStar = pk0[0];
//    rhoStar = 1. / v;
//    vStar = vL;
//    wStar = wL;
//  }
//  //Cas reservoir
//  //-------------
//  else { //Entree reservoir => etat star droite solution
//    //Calcul enthalpie totale etat reservoir
//    double H0(0.);
//    double v0 = 0.;
//    H0 = eos->calculEnthalpieTotale(rhok0[0], pk0[0], v0);
//
//    //PORCEDE ITERATIF CALCUL DE LA PRESSION SOLUTION
//    //-----------------------------------------------
//    int iteration(0);
//    double p(0.5*pk0[0]);
//    double f(0.), df(1.);
//    double dv, h, dh, drho;
//    double uStarR(0.), duStarR(0.), uStarL(0.), duStarL(0.);
//    do {
//      p -= f / df; iteration++;
//      if (iteration > 50) Erreurs::messageErreur("nombre iterations trop grand dans resolRiemannRes modEuler");
//      //Verification pression physique ?
//      eos->verifieEtCorrigePression(p);
//      if (p > pk0[0]) { p = pk0[0] - 1e-6; }
//
//      //Relations reservoir a droite (H=cte et s=cste)
//      rhoStar = eos->calculDensiteIsentrope(pk0[0], rhok0[0], p);
//      h = eos->calculEnthalpieIsentrope(pk0[0], rhok0[0], p, &dh);
//      uStarR = -sqrt(2.*(H0 - h));
//      duStarR = -dh / uStarR; ;
//
//      //Relations a gauche (isentropiques //FP//DEV// chocs a mettre)
//      double dmL;
//      v = 1.0 / eos->calculDensiteIsentrope(pL, rhoL, p, &drho); //FP//DEV// mettre hugoniot
//      //v = 1.0 / eos->calculDensiteHugoniot(pL, rhoL, p, &drho);
//      dv = - v*v*drho;
//      vmv0 = v - 1. / rhoL;
//      if (fabs(vmv0) > 1e-10) {
//        mL = sqrt((pL - p) / vmv0);
//        dmL = 0.5*(-vmv0 + (p - pL)*dv) / (vmv0*vmv0) / mL;
//      }
//      else { //si depassement densite limite sous choc => relations acoustiques
//        mL = zL;
//        dmL = 0.;
//      }
//      sL = uL - mL / rhoL;
//      if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
//      uStarL = uL + mL*vmv0;
//      duStarL = dmL*vmv0 + mL*dv;
//
//      //Fonction a resoudre
//      f = uStarR - uStarL;
//      df = duStarR - duStarL;
//
//    } while (fabs(f)>1e-3);
//
//    pStar = p;
//    uStar = 0.5*(uStarL + uStarR);
//    vStar = 0.;
//    wStar = 0.;
//  }
//
//  eStar = eos->calculEnergie(rhoStar, pStar);
//
//  fluxTempEulerHomogene.m_masse = rhoStar*uStar;
//  fluxTempEulerHomogene.m_qdm.setX(rhoStar*uStar*uStar + pStar);
//  fluxTempEulerHomogene.m_qdm.setY(rhoStar*uStar*vStar);
//  fluxTempEulerHomogene.m_qdm.setZ(rhoStar*uStar*wStar);
//  fluxTempEulerHomogene.m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vStar*vStar + wStar*wStar)) + pStar)*uStar; //FP//DEV// Vitesse transverse nulles dans le reservoir, non ?
//
//  //vitesse bord
//  fluxTempEulerHomogene.m_sM = uStar;
//}
//
////****************************************************************************
//
//void ModEulerHomogene::resolRiemannSortie(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double p0, double *debitSurf) const
//{
//  Eos *eos;
//
//  double cL, sL, zL;
//  double uL, pL, rhoL, vL, wL;
//  double uStar(0.), rhoStar(0.), pStar(0.), eStar(0.);
//
//  double *dataEos;
//
//  Phase *phaseGauche(0);
//  phaseGauche = cellGauche.getPhase(0);
//
//  eos = phaseGauche->getEos();
//  uL = phaseGauche->getU();
//  vL = phaseGauche->getV();
//  wL = phaseGauche->getW();
//  pL = phaseGauche->getPression();
//  rhoL = phaseGauche->getDensite();
//  cL = eos->calculVitesseSon(rhoL, pL);
//  zL = rhoL*cL;
//
//  sL = uL - cL;
//  if (fabs(sL)>1.e-3) dtMax = min(dtMax, dxGauche / fabs(sL));
//
//  eos->renvoiInfo(dataEos); //FP//TODO// A enlever : tres tres moche + probleme de desallocation !
//  pStar = p0;
//  rhoStar = rhoL*pow((pStar / pL), 1. / dataEos[0]);
//  uStar = (pL + zL*uL - pStar) / zL;
//
//  //Cas patologique sL>0            //FP//Q// Peut etre cas u<0 a voir aussi
//  if (sL >= 0.) { //Cas sortie supersonique => etat gauche solution
//    uStar = uL;
//    pStar = pL;
//    rhoStar = rhoL;
//  }
//
//  eStar = eos->calculEnergie(rhoStar, pStar);
//  fluxTempEulerHomogene.m_masse = rhoStar*uStar;
//  fluxTempEulerHomogene.m_qdm.setX(rhoStar*uStar*uStar + pStar);
//  fluxTempEulerHomogene.m_qdm.setY(rhoStar*uStar*vL);
//  fluxTempEulerHomogene.m_qdm.setZ(rhoStar*uStar*wL);
//  fluxTempEulerHomogene.m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vL*vL + wL*wL)) + pStar)*uStar;
//
//  //vitesse bord
//  fluxTempEulerHomogene.m_sM = uStar;
//
//  //Sauvegarde du debit surfacique
//  debitSurf[0] = fluxTempEulerHomogene.m_masse;
//}

//****************************************************************************

double ModEulerHomogene::getSM()
{
  return fluxTempEulerHomogene.m_sM;
}

//****************************************************************************

void ModEulerHomogene::projectionRepereAbsolu(const Coord normale, const Coord tangente, const Coord binormale) const
{
  Coord fluxProjete;
  fluxProjete.setX(normale.getX()*fluxTempEulerHomogene.m_qdm.getX() + tangente.getX()*fluxTempEulerHomogene.m_qdm.getY() + binormale.getX()*fluxTempEulerHomogene.m_qdm.getZ());
  fluxProjete.setY(normale.getY()*fluxTempEulerHomogene.m_qdm.getX() + tangente.getY()*fluxTempEulerHomogene.m_qdm.getY() + binormale.getY()*fluxTempEulerHomogene.m_qdm.getZ());
  fluxProjete.setZ(normale.getZ()*fluxTempEulerHomogene.m_qdm.getX() + tangente.getZ()*fluxTempEulerHomogene.m_qdm.getY() + binormale.getZ()*fluxTempEulerHomogene.m_qdm.getZ());
  fluxTempEulerHomogene.m_qdm.setXYZ(fluxProjete.getX(), fluxProjete.getY(), fluxProjete.getZ());
}

//****************************************************************************

int ModEulerHomogene::getLiq()
{
  return m_liq;
}

//****************************************************************************

int ModEulerHomogene::getVap()
{
  return m_vap;
}

//****************************************************************************

string ModEulerHomogene::quiSuisJe() const
{
  return m_nom;
}

//****************************************************************************