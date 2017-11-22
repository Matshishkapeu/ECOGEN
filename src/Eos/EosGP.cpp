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
#include <cstring>
#include "EosGP.h"

using namespace std;

//***********************************************************************

EosGP::EosGP(){}

//***********************************************************************

EosGP::EosGP(vector<string> &nomParametresEos, int &numero) :
    Eos(numero)
{
  nomParametresEos.push_back("gamma");
  nomParametresEos.push_back("cv");
  nomParametresEos.push_back("energieRef");
  nomParametresEos.push_back("entropieRef");
}

//***********************************************************************

//Constructeur GazParfait ameliore
EosGP::EosGP(string nom, double gamma, double cv, double eRef, int &numero) :
    Eos(nom, numero),
    m_gamma(gamma),
    m_cv(cv),
    m_eRef(eRef),
    m_sRef(0.)
{}

//***********************************************************************

EosGP::~EosGP(){}

//***********************************************************************

//Mod
void EosGP::renvoiInfo(double *&data) const
{
	int nombre = 3;
	data = new double[nombre];

	data[0] = m_gamma;
	data[1] = m_cv;
	data[2] = m_eRef;
}

//***********************************************************************

void EosGP::attributParametresEos(string nom, vector<double> parametresEos)
{
  m_nom   = nom;
  assert(parametresEos.size() == 4);
  m_gamma = parametresEos[0];
  m_cv    = parametresEos[1];
  m_eRef  = parametresEos[2];
  m_sRef  = parametresEos[3];
}

//***********************************************************************

//Methodes constantes
//*******************
double EosGP::calculTemperature(const double &densite, const double &pression) const
{
  return pression/(m_gamma-1.)/densite/m_cv;
}

//***********************************************************************

double EosGP::calculEnergie(const double &densite, const double &pression) const
{
  return pression/(m_gamma-1.)/densite + m_eRef;
}

//***********************************************************************

double EosGP::calculPression(const double &densite, const double &energie) const
{
  return (m_gamma-1.)*densite*(energie-m_eRef);
}

//***********************************************************************

double EosGP::calculVitesseSon(const double &densite, const double &pression) const
{
  return sqrt(m_gamma*pression/densite);
}

//***********************************************************************

double EosGP::calculPressionIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &densiteFinale) const
{
  return pressionInitiale*pow(densiteFinale/densiteInitiale,m_gamma);
}

//***********************************************************************

double EosGP::calculPressionHugoniot(const double &pressionInitiale, const double &densiteInitiale, const double &densiteFinale) const
{
  return pressionInitiale*((m_gamma+1.)*densiteFinale-(m_gamma-1.)*densiteInitiale)/((m_gamma+1.)*densiteInitiale-(m_gamma-1.)*densiteFinale);
}

//***********************************************************************

double EosGP::calculDensiteIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *drhodp) const
{
  double densiteFinale(densiteInitiale*pow(pressionFinale/pressionInitiale,1./m_gamma));
  if(drhodp!=NULL) *drhodp = densiteFinale/(m_gamma*pressionFinale);
  return densiteFinale;
}

//***********************************************************************

double EosGP::calculDensiteHugoniot(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *drhodp) const
{
  double num((m_gamma + 1.)*pressionFinale + (m_gamma - 1.)*pressionInitiale);
  double denom((m_gamma - 1.)*pressionFinale + (m_gamma + 1.)*pressionInitiale);
  double densiteFinale(densiteInitiale*num / denom);
  if (drhodp != NULL) *drhodp = densiteInitiale*2.*(m_gamma + 1.)*pressionInitiale / (denom*denom);
  return densiteFinale;
}

//***********************************************************************

double EosGP::calculEnthalpieIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *dhdp) const
{
  double rhoFinale, drho;
  rhoFinale = this->calculDensiteIsentrope(pressionInitiale, densiteInitiale, pressionFinale, &drho);
  double enthalpieFinale(m_gamma*pressionFinale/(m_gamma-1.)/rhoFinale+m_eRef);
  if (dhdp != NULL) *dhdp = m_gamma / (m_gamma - 1.)*(rhoFinale - pressionFinale*drho) / (rhoFinale*rhoFinale);
  return enthalpieFinale;
}

//***********************************************************************

double EosGP::calculDensiteSaturation(const double &pression, const double &Tsat, const double &dTsatdP, double *drhodp) const
{
  double rho;
  if (drhodp != NULL) {
    *drhodp = (m_gamma - 1.)*m_cv*Tsat - pression*(m_gamma - 1.)*m_cv*dTsatdP;
    *drhodp /= (((m_gamma - 1.)*m_cv*Tsat)*((m_gamma - 1.)*m_cv*Tsat));
  }
  rho = pression/((m_gamma - 1.)*m_cv*Tsat);
  return rho;
}

//***********************************************************************

double EosGP::calculRhoEnergieSaturation(const double &pression, const double &rho, const double &drhodp, double *drhoedp) const
{
  double rhoe;
  if (drhoedp != NULL) {  *drhoedp = 1. / (m_gamma - 1.) + drhodp*m_eRef; }
  rhoe = pression / (m_gamma - 1.) + rho*m_eRef;
  return rhoe;
}


//***********************************************************************

void EosGP::renvoiSpecialEosMelange(double &gamPinfSurGamMoinsUn, double &eRef, double &unSurGamMoinsUn) const
{
  gamPinfSurGamMoinsUn = 0.;
  eRef = m_eRef;
  unSurGamMoinsUn = 1. / (m_gamma - 1.);
}

//***********************************************************************

double EosGP::vfpfh(const double &pression, const double &enthalpie) const
{
  return (m_gamma - 1.)*(enthalpie - m_eRef) / (m_gamma*pression);
}

//***********************************************************************

double EosGP::dvdpch(const double &pression, const double &enthalpie) const
{
  return (1. - m_gamma) / m_gamma * (enthalpie - m_eRef) / (pression*pression);
}

//***********************************************************************

double EosGP::dvdhcp(const double &pression, const double &enthalpie) const
{
  return (m_gamma - 1.) / m_gamma / pression;
}

//***********************************************************************

void EosGP::verifiePression(const double &pression) const
{
  if (pression <= 1.e-6) erreurs.push_back(Erreurs("pression trop faible dans EosGP"));
}

//***********************************************************************

void EosGP::verifieEtCorrigePression(double &pression) const
{
  if (pression <= 1.e-6) pression = 1.e-6;
}

//***********************************************************************

double EosGP::getGamma() const { return m_gamma; }

//***********************************************************************

double EosGP::getCv() const { return m_cv; }

//***********************************************************************

double EosGP::getERef() const { return m_eRef; }

//***********************************************************************

double EosGP::getSRef() const { return m_sRef; }


//***********************************************************************