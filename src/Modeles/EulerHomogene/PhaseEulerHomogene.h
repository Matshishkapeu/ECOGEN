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

#ifndef PHASEEULERHOMOGENE_H
#define PHASEEULERHOMOGENE_H

#include "../Phase.h"
#include "../../Eos/Eos.h"
#include <fstream>

class PhaseEulerHomogene : public Phase
{
public:
  PhaseEulerHomogene();
  PhaseEulerHomogene(tinyxml2::XMLElement *materiau, Eos *eos, std::string nomFichier);
  PhaseEulerHomogene(double alpha, double densite, double pression, double u, double v, double w, Eos *eos);
  virtual ~PhaseEulerHomogene();

  virtual void ecritPhase(std::ofstream &fluxFichier) const;
  virtual void alloueEtCopiePhase(Phase **vecPhase);
  virtual void copiePhase(Phase &vecPhase);
  virtual void calculsEtendusPhases(const Coord &vitesse); //Pour remplir les donnees additionelles

  virtual void projection(const Coord &normale, const Coord &tangente, const Coord &binormale) {};
  virtual void projectionRepereAbsolu(const Coord &normale, const Coord &tangente, const Coord &binormale) {};

  //Accesseurs donnees
  virtual double getAlpha() const;
  virtual double getDensite() const;
  virtual double getPression() const;
  virtual double getU() const { return 0.; };
  virtual double getV() const { return 0.; };
  virtual double getW() const { return 0.; };
  virtual Coord getVitesse() const { return 0; };
  virtual Eos* getEos() const;
  virtual double getEnergie() const;
  virtual double getVitesseSon() const;
  virtual double getEnergieTotale() const;
  virtual double getTemperature() const;

  virtual void setAlpha(double alpha);
  virtual void setDensite(double densite);
  virtual void setPression(double pression);
  virtual void setVitesse(const double &u, const double &v, const double &w) {};
  virtual void setVitesse(const Coord &vit) {};
  virtual void setU(const double &u) {};
  virtual void setV(const double &v) {};
  virtual void setW(const double &w) {};
  virtual void setEos(Eos *eos);
  virtual void setEnergie(double energie);
  virtual void setVitesseSon(double vitesseSon);
  virtual void setEnergieTotale(double energieTotale);

  //Specifique a l ecriture des donnees
  virtual int getNombreScalaires() const { return 3; };
  virtual int getNombreVecteurs() const { return 0; };
  virtual double renvoieScalaire(const int &numVar) const;
  virtual Coord* renvoieVecteur(const int &numVar) { return 0; };
  virtual std::string renvoieNomScalaire(const int &numVar) const;
  virtual std::string renvoieNomVecteur(const int &numVar) const { return 0; };

  //Specifique a la reprise depuis Fichier
  virtual void setScalaire(const int &numVar, const double &valeur);

  //Specifique au parallele
  virtual int nombreVariablesATransmettre() const;
  virtual void rempliTampon(double *tampon, int &compteur) const;
  virtual void recupereTampon(double *tampon, int &compteur, Eos **eos);

  //Specifique a l ordre 2
  virtual void calculPentesPhase(const Phase &pGauche, const Phase &pDroite, const double &distance);
  virtual void miseAZero();
  virtual void extrapole(const Phase &pente, const double &distance);
  virtual void limitePentes(const Phase &penteGauche, const Phase &penteDroite, Limiteur &limiteurGlobal, Limiteur &limiteurInterface);

	//Specifique a l'ordre 2 parallele
	virtual int nombrePentesATransmettre() const;
	virtual void rempliTamponPentes(double *tampon, int &compteur) const;
	virtual void recupereTamponPentes(double *tampon, int &compteur);

  //Verifications
  virtual void verifiePhase(const std::string &message = "") const;
  virtual void verifieEtCorrigePhase(const std::string &message = "");

  //Operateurs
  virtual void changeSigne();
  virtual void multiplieEtAjoute(const Phase &pentesPhasesTemp, const double &coeff);
  virtual void diviser(const double &coeff);

protected:
  double m_alpha;
  double m_densite;
  double m_pression;
  Eos *m_eos;
  double m_energie;
  double m_vitesseSon;
  double m_energieTotale;
private:
};

#endif // PHASEEULERHOMOGENE_H
