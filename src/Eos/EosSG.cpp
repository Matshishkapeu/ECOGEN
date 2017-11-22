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
#include "EosSG.h"

using namespace std;

//***********************************************************************

EosSG::EosSG(){}

//***********************************************************************

EosSG::EosSG(vector<string> &nomParametresEos, int &numero) :
    Eos(numero)
{
  nomParametresEos.push_back("gamma");
  nomParametresEos.push_back("pInf");
  nomParametresEos.push_back("cv");
  nomParametresEos.push_back("energieRef");
  nomParametresEos.push_back("entropieRef");
}

//***********************************************************************

//Constructeur StiffenedGas ameliore
EosSG::EosSG(string nom, double gamma, double pInf, double cv, double eRef, int &numero) :
    Eos(nom, numero),
    m_gamma(gamma),
    m_pInf(pInf),
    m_cv(cv),
    m_eRef(eRef),
    m_sRef(0.)
{}

//***********************************************************************

EosSG::~EosSG(){}

//***********************************************************************

//Mod
void EosSG::renvoiInfo(double *&data) const
{
	int nombre = 4;
	data = new double[nombre];

	data[0] = m_gamma;
	data[1] = m_pInf;
	data[2] = m_cv;
	data[3] = m_eRef;
}

//***********************************************************************

void EosSG::attributParametresEos(string nom, vector<double> parametresEos)
{
  m_nom = nom;
  assert(parametresEos.size() == 5);
  m_gamma = parametresEos[0];
  m_pInf  = parametresEos[1];
  m_cv    = parametresEos[2];
  m_eRef  = parametresEos[3];
  m_sRef  = parametresEos[4];
}

//***********************************************************************

//Methodes constantes
//*******************
double EosSG::calculTemperature(const double &densite, const double &pression) const
{
  return (pression+m_pInf)/(m_gamma-1.)/densite/m_cv;
}

//***********************************************************************

double EosSG::calculEnergie(const double &densite, const double &pression) const
{
  return (pression+m_gamma*m_pInf)/(m_gamma-1.)/densite + m_eRef;
}

//***********************************************************************

double EosSG::calculPression(const double &densite, const double &energie) const
{
  return (m_gamma-1.)*densite*(energie-m_eRef)-m_gamma*m_pInf;
}

//***********************************************************************

double EosSG::calculVitesseSon(const double &densite, const double &pression) const
{
  return sqrt(m_gamma*(pression+m_pInf)/densite);
}

//***********************************************************************

double EosSG::calculPressionIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &densiteFinale) const
{
  return (pressionInitiale+m_pInf)*pow(densiteFinale/densiteInitiale,m_gamma)-m_pInf;
}

//***********************************************************************

double EosSG::calculPressionHugoniot(const double &pressionInitiale, const double &densiteInitiale, const double &densiteFinale) const
{
  return (pressionInitiale+m_pInf)*((m_gamma+1.)*densiteFinale-(m_gamma-1.)*densiteInitiale)/((m_gamma+1.)*densiteInitiale-(m_gamma-1.)*densiteFinale)-m_pInf;
}

//***********************************************************************

double EosSG::calculDensiteIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *drhodp) const
{
  double densiteFinale(densiteInitiale*pow((pressionFinale+m_pInf)/(pressionInitiale+m_pInf),1./m_gamma));
  if (drhodp != NULL) *drhodp = densiteFinale/(m_gamma*(pressionFinale+m_pInf));
  return densiteFinale;
}

//***********************************************************************

double EosSG::calculDensiteHugoniot(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *drhodp) const
{
  double num((m_gamma+1.)*(pressionFinale+m_pInf)+ (m_gamma - 1.)*(pressionInitiale + m_pInf));
  double denom((m_gamma - 1.)*(pressionFinale + m_pInf) + (m_gamma + 1.)*(pressionInitiale + m_pInf));
  double densiteFinale(densiteInitiale*num/denom);
  if (drhodp != NULL) *drhodp = densiteInitiale*2.*(m_gamma+1.)*(pressionInitiale+m_pInf)/(denom*denom);
  return densiteFinale;
}

//***********************************************************************

double EosSG::calculEnthalpieIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *dhdp) const
{
  double rhoFinale, drho;
  rhoFinale = this->calculDensiteIsentrope(pressionInitiale, densiteInitiale, pressionFinale, &drho);
  double enthalpieFinale(m_gamma*(pressionFinale+m_pInf) / (m_gamma - 1.) / rhoFinale + m_eRef);
  if (dhdp != NULL) *dhdp = m_gamma / (m_gamma - 1.)*(rhoFinale - (pressionFinale+m_pInf)*drho) / (rhoFinale*rhoFinale);
  return enthalpieFinale;
}

//***********************************************************************

double EosSG::calculDensiteSaturation(const double &pression, const double &Tsat, const double &dTsatdP, double *drhodp) const
{
  double rho;
  if (drhodp != NULL) {
    *drhodp = (m_gamma - 1.)*m_cv*Tsat - (pression + m_pInf)*(m_gamma - 1.)*m_cv*dTsatdP;
    *drhodp /= (((m_gamma - 1.)*m_cv*Tsat)*((m_gamma - 1.)*m_cv*Tsat));
  }
  rho = (pression + m_pInf)/((m_gamma - 1.)*m_cv*Tsat);
  return rho;
}

//***********************************************************************

double EosSG::calculRhoEnergieSaturation(const double &pression, const double &rho, const double &drhodp, double *drhoedp) const
{
  double rhoe;
  if (drhoedp != NULL) { *drhoedp = 1./(m_gamma-1.)+drhodp*m_eRef; }
  rhoe = (pression + m_gamma*m_pInf) / (m_gamma - 1.) + rho*m_eRef;
  return rhoe;
}

//***********************************************************************

void EosSG::renvoiSpecialEosMelange(double &gamPinfSurGamMoinsUn, double &eRef, double &unSurGamMoinsUn) const
{
  gamPinfSurGamMoinsUn = m_gamma*m_pInf/(m_gamma-1.);
  eRef = m_eRef;
  unSurGamMoinsUn = 1. / (m_gamma - 1.);
}

//***********************************************************************

double EosSG::vfpfh(const double &pression, const double &enthalpie) const
{
  return (m_gamma - 1.)*(enthalpie-m_eRef) / (m_gamma*(pression+m_pInf));
}

//***********************************************************************

double EosSG::dvdpch(const double &pression, const double &enthalpie) const
{
  return (1. - m_gamma) / m_gamma * (enthalpie - m_eRef) / ((pression + m_pInf)*(pression + m_pInf));
}

//***********************************************************************

double EosSG::dvdhcp(const double &pression, const double &enthalpie) const
{
  return (m_gamma - 1.) / m_gamma / (pression + m_pInf);
}

//***********************************************************************

void EosSG::verifiePression(const double &pression) const
{
  if (pression <= -(1. - 1.e-8)*m_pInf + 1.e-8) erreurs.push_back(Erreurs("pression trop faible dans EosSG"));
}

//***********************************************************************

void EosSG::verifieEtCorrigePression(double &pression) const
{
  if (pression <= -(1. - 1.e-8)*m_pInf + 1.e-8) pression = -(1. - 1.e-8)*m_pInf + 1.e-8;
}

//***********************************************************************

double EosSG::getGamma() const { return m_gamma; }

//***********************************************************************

double EosSG::getPInf() const { return m_pInf; }

//***********************************************************************

double EosSG::getCv() const{ return m_cv; }

//***********************************************************************

double EosSG::getERef() const { return m_eRef; }

//***********************************************************************

double EosSG::getSRef() const { return m_sRef; }

//***********************************************************************