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

#ifndef PHASE_H
#define PHASE_H

#include <fstream>
#include "../Erreurs.h"
#include "../Eos/Eos.h"
#include "../Maths/Coord.h"
#include "../libTierces/tinyxml2.h"
#include "../Ordre2/EnteteLimiteur.h"

enum Prim { vecPhases, vecPhasesO2, vecPentes };

class Phase
{
  public:
    Phase();
    virtual ~Phase();

    virtual void ecritPhase(std::ofstream &fluxFichier) const { Erreurs::messageErreur("ecritPhase non accessible pour le type de phase demande"); };
    virtual void alloueEtCopiePhase(Phase **vecPhase) { Erreurs::messageErreur("alloueEtCopiePhase non accessible pour le type de phase demande"); };
    virtual void copiePhase(Phase &vecPhase) { Erreurs::messageErreur("copiePhase non accessible pour le type de phase demande"); };
    virtual void calculsEtendusPhases(const Coord &vitesse) { Erreurs::messageErreur("calculsEtendusPhases non accessible pour le type de phase demande"); };
    
    virtual void projection(const Coord &normale, const Coord &tangente, const Coord &binormale) { Erreurs::messageErreur("projection non accessible pour le type de phase demande"); };
    virtual void projectionRepereAbsolu(const Coord &normale, const Coord &tangente, const Coord &binormale) { Erreurs::messageErreur("projectionRepereAbsolu non accessible pour le type de phase demande"); };

    //Accesseurs donnees
    virtual double getAlpha() const { Erreurs::messageErreur("getAlpha non accessible pour le type de phase demande"); return 0.; };
    virtual double getDensite() const { Erreurs::messageErreur("getDensite non accessible pour le type de phase demande"); return 0.; };
    virtual double getPression() const { Erreurs::messageErreur("getPression non accessible pour le type de phase demande"); return 0.; };
    virtual double getU() const { Erreurs::messageErreur("getU non accessible pour le type de phase demande"); return 0.; };
    virtual double getV() const { Erreurs::messageErreur("getV non accessible pour le type de phase demande"); return 0.; };
    virtual double getW() const { Erreurs::messageErreur("getW non accessible pour le type de phase demande"); return 0.; };
    virtual Coord getVitesse() const { Erreurs::messageErreur("getVitesse non accessible pour le type de phase demande"); return 0; };
    virtual Eos* getEos() const { Erreurs::messageErreur("EOS non accessible pour le type de phase demande"); return 0; };
    virtual double getEnergie() const { Erreurs::messageErreur("getEnergie impossible avec type de phase demande"); return 0.; };
    virtual double getVitesseSon() const { Erreurs::messageErreur("getVitesseSon impossible avec type de phase demande"); return 0.; };
    virtual double getEnergieTotale() const { Erreurs::messageErreur("getEnergieTotale impossible avec type de phase demande"); return 0.; };
    virtual double getTemperature() const { Erreurs::messageErreur("getT impossible avec type de phase demande"); return 0.; };

    virtual void setAlpha(double alpha) { Erreurs::messageErreur("setAlpha non implemente pour phase utilise"); };
    virtual void setDensite(double densite) { Erreurs::messageErreur("setDensite non implemente pour phase utilise"); };
    virtual void setPression(double pression) { Erreurs::messageErreur("setPression non implemente pour phase utilise"); };
    virtual void setVitesse(const double &u, const double &v, const double &w) { Erreurs::messageErreur("setVitesse non implemente pour phase utilise"); };
    virtual void setVitesse(const Coord &vit) { Erreurs::messageErreur("setVitesse non implemente pour phase utilise"); };
    virtual void setU(const double &u) { Erreurs::messageErreur("setU non implemente pour phase utilise"); };
    virtual void setV(const double &v) { Erreurs::messageErreur("setV non implemente pour phase utilise"); };
    virtual void setW(const double &w) { Erreurs::messageErreur("setW non implemente pour phase utilise"); };
    virtual void setEos(Eos *eos) { Erreurs::messageErreur("Impossible d associer une EOS au type de phase demande"); };
    virtual void setEnergie(double energie) { Erreurs::messageErreur("setEnergie non implemente pour phase utilise"); };
    virtual void setVitesseSon(double vitesseSon) { Erreurs::messageErreur("setVitesseSon non implemente pour phase utilise"); };
    virtual void setEnergieTotale(double energieTotale) { Erreurs::messageErreur("setEnergieTotale non implemente pour phase utilise"); };

    //Specifique a l ecriture des donnees
    virtual int getNombreScalaires() const { Erreurs::messageErreur("getNombreScalaires non implemente pour phase utilise"); return 0; };
    virtual int getNombreVecteurs() const { Erreurs::messageErreur("getNombreVecteurs non implemente pour phase utilise"); return 0; };
    virtual double renvoieScalaire(const int &numVar) const { Erreurs::messageErreur("renvoieScalaire non implemente pour phase utilise"); return 0.; };
    virtual Coord* renvoieVecteur(const int &numVar) { Erreurs::messageErreur("renvoieVecteur non implemente pour phase utilise"); return 0; };
    virtual std::string renvoieNomScalaire(const int &numVar) const { Erreurs::messageErreur("renvoieNomScalaire non implemente pour phase utilise"); return 0; };
    virtual std::string renvoieNomVecteur(const int &numVar) const { Erreurs::messageErreur("renvoieNomVecteur non implemente pour phase utilise"); return 0; };

    //Specifique a la reprise depuis Fichier
    virtual void setScalaire(const int &numVar, const double &valeur) { Erreurs::messageErreur("setScalaire non implemente pour phase utilise"); };

    //Specifique au parallele
    virtual int nombreVariablesATransmettre() const { Erreurs::messageErreur("nombreVariablesATransmettre non implemente pour phase utilise"); return 0; };
    virtual void rempliTampon(double *tampon, int &compteur) const { Erreurs::messageErreur("rempliTampon non implemente pour phase utilise"); };
    virtual void recupereTampon(double *tampon, int &compteur, Eos **eos) { Erreurs::messageErreur("recupereTampon non implemente pour phase utilise"); };

    //Specifique a l ordre 2
    virtual void calculPentesPhase(const Phase &pGauche, const Phase &pDroite, const double &distance) { Erreurs::messageErreur("calculPentesPhase non implemente pour phase utilise"); };
    virtual void miseAZero() { Erreurs::messageErreur("miseAZero non implemente pour phase utilise"); };
    virtual void extrapole(const Phase &pente, const double &distance) { Erreurs::messageErreur("extrapole non implemente pour phase utilise"); };
    virtual void limitePentes(const Phase &penteGauche, const Phase &penteDroite, Limiteur &limiteurGlobal, Limiteur &limiteurInterface) { Erreurs::messageErreur("limitePentes non implemente pour phase utilise"); };

		//Specifique a l ordre 2 parallele
		virtual int nombrePentesATransmettre() const { Erreurs::messageErreur("nombrePentesATransmettre non implemente pour phase utilise"); return 0; };
		virtual void rempliTamponPentes(double *tampon, int &compteur) const { Erreurs::messageErreur("rempliTamponPentes non implemente pour phase utilise"); };
		virtual void recupereTamponPentes(double *tampon, int &compteur) { Erreurs::messageErreur("recupereTamponPentes non implemente pour phase utilise"); };

    //Verifications
    virtual void verifiePhase(const std::string &message = "") const { Erreurs::messageErreur("verifiePhase non prevu pour type de phase"); };
    virtual void verifieEtCorrigePhase(const std::string &message = "") { Erreurs::messageErreur("verifieEtCorrigePhase non prevu pour type de phase"); };

    //Operateurs
    virtual void changeSigne() { Erreurs::messageErreur("changeSigne non implemente pour phase utilise"); };
    virtual void multiplieEtAjoute(const Phase &pentesPhasesTemp, const double &coeff) { Erreurs::messageErreur("multiplieEtAjoute non implemente pour phase utilise"); };
    virtual void diviser(const double &coeff) { Erreurs::messageErreur("diviser non implemente pour phase utilise"); };

  protected:

};

#endif // PHASE_H
