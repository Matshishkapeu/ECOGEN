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

#include "PhaseEulerHomogene.h"
#include "../../Eos/Eos.h"
#include <fstream>

using namespace std;
using namespace tinyxml2;

//***************************************************************************

PhaseEulerHomogene::PhaseEulerHomogene() :m_alpha(1.0), m_densite(0.), m_pression(0.), m_eos(0), m_energie(0.), m_energieTotale(0.), m_vitesseSon(0.) {}

//***************************************************************************

/*! Constructeur materiau a partir d une lecture au format XML
! ex : 			<donneesFluide alpha="1.0" densite="1.0" pression="1e5">
!              <vitesse x = "0." y = "0." z = "0." />
!           </donneesFluide>
*/
PhaseEulerHomogene::PhaseEulerHomogene(XMLElement *materiau, Eos *eos, string nomFichier) : m_eos(eos), m_energie(0.), m_energieTotale(0.), m_vitesseSon(0.)
{
  XMLElement *sousElement(materiau->FirstChildElement("donneesFluide"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesFluide", nomFichier, __FILE__, __LINE__);
  //Recuperation des donnes
  //-----------------------
  XMLError erreur;
  //alpha
  erreur = sousElement->QueryDoubleAttribute("alpha", &m_alpha);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("alpha", nomFichier, __FILE__, __LINE__);
  //densite
  erreur = sousElement->QueryDoubleAttribute("densite", &m_densite);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("densite", nomFichier, __FILE__, __LINE__);
  //pression
  erreur = sousElement->QueryDoubleAttribute("pression", &m_pression);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("pression", nomFichier, __FILE__, __LINE__);
}

//***************************************************************************

PhaseEulerHomogene::PhaseEulerHomogene(double alpha, double densite, double pression, double u, double v, double w, Eos *eos) :
  m_alpha(alpha), m_densite(densite), m_pression(pression), m_eos(eos) {}

//***************************************************************************

PhaseEulerHomogene::~PhaseEulerHomogene() {}

//***************************************************************************

void PhaseEulerHomogene::ecritPhase(ofstream &fluxFichier) const
{
  fluxFichier << m_alpha << " " << m_densite << " " << m_pression << " ";
  //fluxFichier << m_alpha << " " << m_densite << " " << m_pression << " " << this->getTemperature() << " ";
}

//***************************************************************************

void PhaseEulerHomogene::alloueEtCopiePhase(Phase **vecPhase)
{
  *vecPhase = new PhaseEulerHomogene(*this);
}

//***************************************************************************

void PhaseEulerHomogene::copiePhase(Phase &phase)
{
  m_alpha = phase.getAlpha();
  m_densite = phase.getDensite();
  m_pression = phase.getPression();
  m_eos = phase.getEos();
  m_energie = phase.getEnergie();
  m_vitesseSon = phase.getVitesseSon();
  m_energieTotale = phase.getEnergieTotale();
}

//***************************************************************************

void PhaseEulerHomogene::calculsEtendusPhases(const Coord &vitesse)
{
  m_energie = m_eos->calculEnergie(m_densite, m_pression);
  m_vitesseSon = m_eos->calculVitesseSon(m_densite, m_pression);
  m_energieTotale = m_energie + 0.5*vitesse.normeCarre();
}

//****************************************************************************
//************************** ACCESSEURS DONNEES ******************************
//****************************************************************************

double PhaseEulerHomogene::getAlpha() const { return m_alpha; }

//***************************************************************************

double PhaseEulerHomogene::getDensite() const { return m_densite; }

//***************************************************************************

double PhaseEulerHomogene::getPression() const { return m_pression; }

//***************************************************************************

Eos* PhaseEulerHomogene::getEos() const { return m_eos; }

//***************************************************************************

double PhaseEulerHomogene::getEnergie() const { return m_energie; }

//***************************************************************************

double PhaseEulerHomogene::getVitesseSon() const { return m_vitesseSon; }

//***************************************************************************

double PhaseEulerHomogene::getEnergieTotale() const { return m_energieTotale; }

//***************************************************************************

double PhaseEulerHomogene::getTemperature() const { return m_eos->calculTemperature(m_densite, m_pression); }

//***************************************************************************

void PhaseEulerHomogene::setAlpha(double alpha) { m_alpha = alpha; }

//***************************************************************************

void PhaseEulerHomogene::setDensite(double densite) { m_densite = densite; }

//***************************************************************************

void PhaseEulerHomogene::setPression(double pression) { m_pression = pression; }

//***************************************************************************

void PhaseEulerHomogene::setEos(Eos *eos) { m_eos = eos; }

//***************************************************************************

void PhaseEulerHomogene::setEnergie(double energie) { m_energie = energie; }

//***************************************************************************

void PhaseEulerHomogene::setVitesseSon(double vitesseSon) { m_vitesseSon = vitesseSon; }

//***************************************************************************

void PhaseEulerHomogene::setEnergieTotale(double energieTotale) { m_energieTotale = energieTotale; }

//****************************************************************************
//************************* ECRITURE DES DONNEES *****************************
//****************************************************************************

double PhaseEulerHomogene::renvoieScalaire(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_alpha; break;
  case 2:
    return m_densite; break;
  case 3:
    return m_pression; break;
  default:
    return 0.; break;
  }
}

//***************************************************************************

string PhaseEulerHomogene::renvoieNomScalaire(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Alpha"; break;
  case 2:
    return "MasseVolumique"; break;
  case 3:
    return "Pression"; break;
  default:
    return "SansNom"; break;
  }
}

//****************************************************************************
//*************************** REPRISE FICHIER ********************************
//****************************************************************************

void PhaseEulerHomogene::setScalaire(const int &numVar, const double &valeur)
{
  switch (numVar)
  {
  case 1:
    m_alpha = valeur; break;
  case 2:
    m_densite = valeur; break;
  case 3:
    m_pression = valeur; break;
  default:
    Erreurs::messageErreur("numVar inconnu dans Phase::setScalaire"); break;
  }
}

//****************************************************************************
//****************************** PARALLELE ***********************************
//****************************************************************************

int PhaseEulerHomogene::nombreVariablesATransmettre() const
{
  //3 variables + numero EOS
  return 4;
}

//***************************************************************************

void PhaseEulerHomogene::rempliTampon(double *tampon, int &compteur) const
{
  tampon[++compteur] = m_alpha;
  tampon[++compteur] = m_densite;
  tampon[++compteur] = m_pression;
  tampon[++compteur] = static_cast<double>(m_eos->getNumero());
}

//***************************************************************************

void PhaseEulerHomogene::recupereTampon(double *tampon, int &compteur, Eos **eos)
{
  m_alpha = tampon[++compteur];
  m_densite = tampon[++compteur];
  m_pression = tampon[++compteur];
  m_eos = eos[static_cast<int>(tampon[++compteur])];
}

//****************************************************************************
//******************************* ORDRE 2 ************************************
//****************************************************************************

void PhaseEulerHomogene::calculPentesPhase(const Phase &pGauche, const Phase &pDroite, const double &distance)
{
  m_alpha = (pDroite.getAlpha() - pGauche.getAlpha()) / distance;
  m_densite = (pDroite.getDensite() - pGauche.getDensite()) / distance;
  m_pression = (pDroite.getPression() - pGauche.getPression()) / distance;
}

//***************************************************************************

void PhaseEulerHomogene::miseAZero()
{
  m_alpha = 0.; m_densite = 0.; m_pression = 0.;
}

//***************************************************************************

void PhaseEulerHomogene::extrapole(const Phase &pente, const double &distance)
{
  m_alpha += pente.getAlpha() * distance;
  m_densite += pente.getDensite() * distance;
  m_pression += pente.getPression() * distance;
}

//***************************************************************************

void PhaseEulerHomogene::limitePentes(const Phase &penteGauche, const Phase &penteDroite, Limiteur &limiteurGlobal, Limiteur &limiteurInterface)
{
  m_alpha = limiteurInterface.limitePente(penteGauche.getAlpha(), penteDroite.getAlpha());
  m_densite = limiteurGlobal.limitePente(penteGauche.getDensite(), penteDroite.getDensite());
  m_pression = limiteurGlobal.limitePente(penteGauche.getPression(), penteDroite.getPression());
}

//****************************************************************************
//************************** ORDRE 2 PARALLELE *******************************
//****************************************************************************

int PhaseEulerHomogene::nombrePentesATransmettre() const
{
	return 3;
}

//***************************************************************************

void PhaseEulerHomogene::rempliTamponPentes(double *tampon, int &compteur) const
{
	tampon[++compteur] = m_alpha;
	tampon[++compteur] = m_densite;
	tampon[++compteur] = m_pression;
}

//***************************************************************************

void PhaseEulerHomogene::recupereTamponPentes(double *tampon, int &compteur)
{
	m_alpha = tampon[++compteur];
	m_densite = tampon[++compteur];
	m_pression = tampon[++compteur];
}

//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseEulerHomogene::verifiePhase(const string &message) const
{
  if (m_alpha <= 1e-10) erreurs.push_back(Erreurs(message + "probleme alpha trop petit dans verifiePhase"));
  if (m_densite <= 1.e-10) erreurs.push_back(Erreurs(message + "probleme masse volumique trop petite dans verifiePhase"));
  m_eos->verifiePression(m_pression);
}

//***************************************************************************

void PhaseEulerHomogene::verifieEtCorrigePhase(const string &message)
{
  if (m_alpha < 1e-10) this->setAlpha(1e-10);
  if (m_alpha > 1. - 1e-10) this->setAlpha(1. - 1e-10);
  if (m_densite < 1.e-10) this->setDensite(1.e-10);
  m_eos->verifieEtCorrigePression(m_pression);
}

//****************************************************************************
//***************************** OPERATEURS ***********************************
//****************************************************************************

void PhaseEulerHomogene::changeSigne()
{
  m_alpha = -m_alpha;
  m_densite = -m_densite;
  m_pression = -m_pression;
}

//***************************************************************************

void PhaseEulerHomogene::multiplieEtAjoute(const Phase &pentesPhasesTemp, const double &coeff)
{
  m_alpha += pentesPhasesTemp.getAlpha()*coeff;
  m_densite += pentesPhasesTemp.getDensite()*coeff;
  m_pression += pentesPhasesTemp.getPression()*coeff;
}

//***************************************************************************

void PhaseEulerHomogene::diviser(const double &coeff)
{
  m_alpha /= coeff;
  m_densite /= coeff;
  m_pression /= coeff;
}
