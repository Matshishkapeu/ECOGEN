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

#include "PhaseEuler.h"
#include "../../Eos/Eos.h"
#include <fstream>

using namespace std;
using namespace tinyxml2;

//***************************************************************************

PhaseEuler::PhaseEuler() :m_densite(0.), m_pression(0.), m_eos(0), m_energie(0.), m_energieTotale(0.), m_vitesseSon(0.)
{
  m_vitesse.setXYZ(0., 0., 0.);
}

//***************************************************************************

/*! Constructeur materiau a partir d une lecture au format XML
  ! ex : 			<donneesFluide densite="1.0" pression="1e5">
  !              <vitesse x = "0." y = "0." z = "0." />
  !           </donneesFluide> 
*/
PhaseEuler::PhaseEuler(XMLElement *materiau, Eos *eos, string nomFichier) : m_eos(eos), m_energie(0.), m_energieTotale(0.), m_vitesseSon(0.)
{
  XMLElement *sousElement(materiau->FirstChildElement("donneesFluide"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesFluide", nomFichier, __FILE__, __LINE__);
  //Recuperation des donnes
  //-----------------------
  XMLError erreur;
  //densite
  erreur = sousElement->QueryDoubleAttribute("densite", &m_densite);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("densite", nomFichier, __FILE__, __LINE__);
  //pression
  erreur = sousElement->QueryDoubleAttribute("pression", &m_pression);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("pression", nomFichier, __FILE__, __LINE__);
  //vitesse
  XMLElement *vitesse(sousElement->FirstChildElement("vitesse"));
  if (vitesse == NULL) throw ErreurXMLElement("vitesse", nomFichier, __FILE__, __LINE__);
  double vitesseX(0.), vitesseY(0.), vitesseZ(0.);
  erreur = vitesse->QueryDoubleAttribute("x", &vitesseX);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("x", nomFichier, __FILE__, __LINE__);
  erreur = vitesse->QueryDoubleAttribute("y", &vitesseY);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("y", nomFichier, __FILE__, __LINE__);
  erreur = vitesse->QueryDoubleAttribute("z", &vitesseZ);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("z", nomFichier, __FILE__, __LINE__);
  m_vitesse.setXYZ(vitesseX, vitesseY, vitesseZ);
}

//***************************************************************************

PhaseEuler::~PhaseEuler(){}

//***************************************************************************

void PhaseEuler::ecritPhase(ofstream &fluxFichier) const
{
  fluxFichier << m_densite << " " << m_pression << " " << m_vitesse.getX() << " ";
  //fluxFichier << m_alpha << " " << m_densite << " " << m_pression << " " << m_vitesse.getX() << " " << m_vitesse.getY() << " " << m_vitesse.getZ() << " " << this->getTemperature() << " ";
}

//***************************************************************************

void PhaseEuler::alloueEtCopiePhase(Phase **vecPhase)
{
  *vecPhase = new PhaseEuler(*this);
}

//***************************************************************************

void PhaseEuler::copiePhase(Phase &phase)
{
  m_densite = phase.getDensite();
  m_pression = phase.getPression();
  m_vitesse = phase.getVitesse();
  m_eos = phase.getEos();
  m_energie = phase.getEnergie();
  m_vitesseSon = phase.getVitesseSon();
  m_energieTotale = phase.getEnergieTotale();
}

//***************************************************************************

void PhaseEuler::calculsEtendusPhases(const Coord &vitesse)
{
  m_energie = m_eos->calculEnergie(m_densite, m_pression);
  m_vitesseSon = m_eos->calculVitesseSon(m_densite, m_pression);
  m_energieTotale = m_energie + 0.5*m_vitesse.normeCarre();
}

//***************************************************************************

void PhaseEuler::projection(const Coord &normale, const Coord &tangente, const Coord &binormale)
{
  m_vitesse.projection(normale, tangente, binormale);
}

//***************************************************************************

void PhaseEuler::projectionRepereAbsolu(const Coord &normale, const Coord &tangente, const Coord &binormale)
{
  m_vitesse.projectionRepereAbsolu(normale, tangente, binormale);
}

//****************************************************************************
//************************** ACCESSEURS DONNEES ******************************
//****************************************************************************

double PhaseEuler::getDensite() const { return m_densite; }

//***************************************************************************

double PhaseEuler::getPression() const { return m_pression; }

//***************************************************************************

double PhaseEuler::getU() const { return m_vitesse.getX(); }
double PhaseEuler::getV() const { return m_vitesse.getY(); }
double PhaseEuler::getW() const { return m_vitesse.getZ(); }

//***************************************************************************

Coord PhaseEuler::getVitesse() const { return m_vitesse; }

//***************************************************************************

Eos* PhaseEuler::getEos() const { return m_eos; }

//***************************************************************************

double PhaseEuler::getEnergie() const { return m_energie; }

//***************************************************************************

double PhaseEuler::getVitesseSon() const { return m_vitesseSon; }

//***************************************************************************

double PhaseEuler::getEnergieTotale() const { return m_energieTotale; }

//***************************************************************************

double PhaseEuler::getTemperature() const { return m_eos->calculTemperature(m_densite, m_pression); }

//***************************************************************************

void PhaseEuler::setDensite(double densite) { m_densite = densite; }

//***************************************************************************

void PhaseEuler::setPression(double pression) { m_pression = pression; }

//***************************************************************************

void PhaseEuler::setVitesse(const double &u, const double &v, const double &w) { m_vitesse.setXYZ(u, v, w); }

//***************************************************************************

void PhaseEuler::setVitesse(const Coord &vit) { m_vitesse = vit; }

//***************************************************************************

void PhaseEuler::setU(const double &u) { m_vitesse.setX(u); }

//***************************************************************************

void PhaseEuler::setV(const double &v) { m_vitesse.setY(v); }

//***************************************************************************

void PhaseEuler::setW(const double &w) { m_vitesse.setZ(w); }

//***************************************************************************

void PhaseEuler::setEos(Eos *eos) { m_eos = eos; }

//***************************************************************************

void PhaseEuler::setEnergie(double energie) { m_energie = energie; }

//***************************************************************************

void PhaseEuler::setVitesseSon(double vitesseSon) { m_vitesseSon = vitesseSon; }

//***************************************************************************

void PhaseEuler::setEnergieTotale(double energieTotale) { m_energieTotale = energieTotale; }

//****************************************************************************
//************************* ECRITURE DES DONNEES *****************************
//****************************************************************************

double PhaseEuler::renvoieScalaire(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_densite; break;
  case 2:
    return m_pression; break;
  default:
    return 0.; break;
  }
}

//***************************************************************************

Coord* PhaseEuler::renvoieVecteur(const int &numVar)
{
  switch (numVar)
  {
  case 1:
    return &m_vitesse; break;
  default:
    return 0; break;
  }
}

//***************************************************************************

string PhaseEuler::renvoieNomScalaire(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return "MasseVolumique"; break;
  case 2:
    return "Pression"; break;
  default:
    return "SansNom"; break;
  }
}

//***************************************************************************

string PhaseEuler::renvoieNomVecteur(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Vitesse"; break;
  default:
    return "SansNom"; break;
  }
}

//****************************************************************************
//*************************** REPRISE FICHIER ********************************
//****************************************************************************

void PhaseEuler::setScalaire(const int &numVar, const double &valeur)
{
  switch (numVar)
  {
  case 1:
    m_densite = valeur; break;
  case 2:
    m_pression = valeur; break;
  default:
    Erreurs::messageErreur("numVar inconnu dans Phase::setScalaire"); break;
  }
}

//****************************************************************************
//****************************** PARALLELE ***********************************
//****************************************************************************

int PhaseEuler::nombreVariablesATransmettre() const
{
  //5 variables + numero EOS
  return 6;
}

//***************************************************************************

void PhaseEuler::rempliTampon(double *tampon, int &compteur) const
{
  tampon[++compteur] = m_densite;
  tampon[++compteur] = m_vitesse.getX();
  tampon[++compteur] = m_vitesse.getY();
  tampon[++compteur] = m_vitesse.getZ();
  tampon[++compteur] = m_pression;
  tampon[++compteur] = static_cast<double>(m_eos->getNumero());
}

//***************************************************************************

void PhaseEuler::recupereTampon(double *tampon, int &compteur, Eos **eos)
{
  m_densite = tampon[++compteur];
  m_vitesse.setX(tampon[++compteur]);
  m_vitesse.setY(tampon[++compteur]);
  m_vitesse.setZ(tampon[++compteur]);
  m_pression = tampon[++compteur];
  m_eos = eos[static_cast<int>(tampon[++compteur])];
}

//****************************************************************************
//******************************* ORDRE 2 ************************************
//****************************************************************************

void PhaseEuler::calculPentesPhase(const Phase &pGauche, const Phase &pDroite, const double &distance)
{
  m_densite = (pDroite.getDensite() - pGauche.getDensite()) / distance;
  m_pression = (pDroite.getPression() - pGauche.getPression()) / distance;
  m_vitesse.setX((pDroite.getVitesse().getX() - pGauche.getVitesse().getX()) / distance);
  m_vitesse.setY((pDroite.getVitesse().getY() - pGauche.getVitesse().getY()) / distance);
  m_vitesse.setZ((pDroite.getVitesse().getZ() - pGauche.getVitesse().getZ()) / distance);
}

//***************************************************************************

void PhaseEuler::miseAZero()
{
  m_densite = 0.; m_pression = 0.;
  m_vitesse.setX(0.); m_vitesse.setY(0.); m_vitesse.setZ(0.);
}

//***************************************************************************

void PhaseEuler::extrapole(const Phase &pente, const double &distance)
{
  m_densite += pente.getDensite() * distance;
  m_pression += pente.getPression() * distance;
  m_vitesse.setX(m_vitesse.getX() + pente.getVitesse().getX() * distance);
  m_vitesse.setY(m_vitesse.getY() + pente.getVitesse().getY() * distance);
  m_vitesse.setZ(m_vitesse.getZ() + pente.getVitesse().getZ() * distance);
}

//***************************************************************************

void PhaseEuler::limitePentes(const Phase &penteGauche, const Phase &penteDroite, Limiteur &limiteurGlobal, Limiteur &limiteurInterface)
{
  m_densite = limiteurGlobal.limitePente(penteGauche.getDensite(), penteDroite.getDensite());
  m_pression = limiteurGlobal.limitePente(penteGauche.getPression(), penteDroite.getPression());
  m_vitesse.setX(limiteurGlobal.limitePente(penteGauche.getVitesse().getX(), penteDroite.getVitesse().getX()));
  m_vitesse.setY(limiteurGlobal.limitePente(penteGauche.getVitesse().getY(), penteDroite.getVitesse().getY()));
  m_vitesse.setZ(limiteurGlobal.limitePente(penteGauche.getVitesse().getZ(), penteDroite.getVitesse().getZ()));
}

//****************************************************************************
//************************** ORDRE 2 PARALLELE *******************************
//****************************************************************************

int PhaseEuler::nombrePentesATransmettre() const
{
	return 5;
}

//***************************************************************************

void PhaseEuler::rempliTamponPentes(double *tampon, int &compteur) const
{
	tampon[++compteur] = m_densite;
	tampon[++compteur] = m_vitesse.getX();
	tampon[++compteur] = m_vitesse.getY();
	tampon[++compteur] = m_vitesse.getZ();
	tampon[++compteur] = m_pression;
}

//***************************************************************************

void PhaseEuler::recupereTamponPentes(double *tampon, int &compteur)
{
	m_densite = tampon[++compteur];
	m_vitesse.setX(tampon[++compteur]);
	m_vitesse.setY(tampon[++compteur]);
	m_vitesse.setZ(tampon[++compteur]);
	m_pression = tampon[++compteur];
}


//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseEuler::verifiePhase(const string &message) const
{
  if (m_densite <= 1.e-10) erreurs.push_back(Erreurs(message + "probleme masse volumique trop petite dans verifiePhase"));
  m_eos->verifiePression(m_pression);
}

//***************************************************************************

void PhaseEuler::verifieEtCorrigePhase(const string &message)
{
  if (m_densite < 1.e-10) this->setDensite(1.e-10);
  m_eos->verifieEtCorrigePression(m_pression);
}

//****************************************************************************
//***************************** OPERATEURS ***********************************
//****************************************************************************

void PhaseEuler::changeSigne()
{
  m_densite = -m_densite;
  m_pression = -m_pression;
  m_vitesse = m_vitesse*-1.;
}

//***************************************************************************

void PhaseEuler::multiplieEtAjoute(const Phase &pentesPhasesTemp, const double &coeff)
{
  m_densite += pentesPhasesTemp.getDensite()*coeff;
  m_pression += pentesPhasesTemp.getPression()*coeff;
  m_vitesse += pentesPhasesTemp.getVitesse()*coeff;
}

//***************************************************************************

void PhaseEuler::diviser(const double &coeff)
{
  m_densite /= coeff;
  m_pression /= coeff;
  m_vitesse /= coeff;
}
