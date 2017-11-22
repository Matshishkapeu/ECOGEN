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

#ifndef MELANGE_H
#define MELANGE_H

#include <vector>

class Melange;

#include "../PhysiqueAdditionnelle/GrandeursPhysAdd.h"

class Melange
{
    public:
      Melange();
      virtual ~Melange();

      virtual void ecritMelange(std::ofstream &fluxFichier) const { Erreurs::messageErreur("ecritMelange non accessible pour le type de melange demande"); };
      virtual void alloueEtCopieMelange(Melange **melange) { Erreurs::messageErreur("alloueEtCopieMelange non accessible pour le type de melange demande"); };
      virtual void alloueYk(const int &nombrePhases) { Erreurs::messageErreur("alloueYk non accessible pour le type de melange demande"); };
      virtual void copieMelange(Melange &melange) { Erreurs::messageErreur("copieMelange non accessible pour le type de melange demande"); };
      virtual double calculDensite(const double *alphak, const double *rhok, const int &nombrePhases) { Erreurs::messageErreur("calculDensite non accessible pour le type de melange demande"); return 0.; };
      virtual double calculPression(const double *alphak, const double *pk, const int &nombrePhases) { Erreurs::messageErreur("calculPression non accessible pour le type de melange demande"); return 0.; };
      virtual double calculEnergieInterne(const double *Yk, const double *ek, const int &nombrePhases) { Erreurs::messageErreur("calculEnergieInterne non accessible pour le type de melange demande"); return 0.; };
      virtual double calculVitesseSonFigee(const double *Yk, const double *ck, const int &nombrePhases) { Erreurs::messageErreur("calculVitesseSonFigee non accessible pour le type de melange demande"); return 0.; };
      
      virtual double calculTsat(const Eos *eosLiq, const Eos *eosVap, const double &pression, double *dTsat=0) { Erreurs::messageErreur("calculTsat non accessible pour le type de melange demande"); return 0.; };

      virtual void calculGrandeursMelange(Phase **vecPhase, const int &nombrePhases) { Erreurs::messageErreur("calculGrandeursMelange non accessible pour le type de melange demande"); };
      virtual void energieInterneVersEnergieTotale(std::vector<GrandeursPhysAdd*> &vecGPA) { Erreurs::messageErreur("energieInterneVersEnergieTotale non accessible pour le type de melange demande"); };
      virtual void energieTotaleVersEnergieInterne(std::vector<GrandeursPhysAdd*> &vecGPA) { Erreurs::messageErreur("energieTotaleVersEnergieInterne non accessible pour le type de melange demande"); };

      virtual void calculEnergieCapillaire() { Erreurs::messageErreur("calculEnergieCapillaire non accessible pour le type de melange demande"); };
      virtual void calculEnergieTotale() { Erreurs::messageErreur("calculEnergieTotale non accessible pour le type de melange demande"); };
      
      virtual void projection(const Coord &normale, const Coord &tangente, const Coord &binormale) { Erreurs::messageErreur("projection non accessible pour le type de melange demande"); };
      virtual void projectionRepereAbsolu(const Coord &normale, const Coord &tangente, const Coord &binormale) { Erreurs::messageErreur("projectionRepereAbsolu non accessible pour le type de melange demande"); };

      //Specifique a l ecriture des donnees
      virtual int getNombreScalaires() const { Erreurs::messageErreur("getNombreScalaires non implemente pour melange utilise"); return 0; };
      virtual int getNombreVecteurs() const { Erreurs::messageErreur("getNombreVecteurs non implemente pour melange utilise"); return 0; };
      virtual double renvoieScalaire(const int &numVar) const { Erreurs::messageErreur("renvoieScalaire non implemente pour melange utilise"); return 0.; };
      virtual Coord* renvoieVecteur(const int &numVar) { Erreurs::messageErreur("renvoieVecteur non implemente pour melange utilise"); return 0; };
      virtual std::string renvoieNomScalaire(const int &numVar) const { Erreurs::messageErreur("renvoieNomScalaire non implemente pour melange utilise"); return 0; };
      virtual std::string renvoieNomVecteur(const int &numVar) const { Erreurs::messageErreur("renvoieNomVecteur non implemente pour melange utilise"); return 0; };

      //Specifique au parallele
      virtual int nombreVariablesATransmettre() const { Erreurs::messageErreur("nombreVariablesATransmettre non implemente pour melange utilise"); return 0; };
      virtual void rempliTampon(double *tampon, int &compteur) const { Erreurs::messageErreur("rempliTampon non implemente pour melange utilise"); };
      virtual void recupereTampon(double *tampon, int &compteur) { Erreurs::messageErreur("recupereTampon non implemente pour melange utilise"); };

      //Specifique a l ordre 2
      virtual void calculPentesMelange(const Melange &pGauche, const Melange &pDroite, const double &distance) { Erreurs::messageErreur("calculPentesMelange non implemente pour melange utilise"); };
      virtual void miseAZero() { Erreurs::messageErreur("miseAZero non implemente pour melange utilise"); };
      virtual void extrapole(const Melange &pente, const double &distance) { Erreurs::messageErreur("extrapole non implemente pour melange utilise"); };
      virtual void limitePentes(const Melange &penteGauche, const Melange &penteDroite, Limiteur &limiteurGlobal) { Erreurs::messageErreur("limitePentes non implemente pour melange utilise"); };

			//Specifique a l'ordre 2 parallele
			virtual int nombrePentesATransmettre() const { Erreurs::messageErreur("nombrePentesATransmettre non implemente pour melange utilise"); return 0; };
			virtual void rempliTamponPentes(double *tampon, int &compteur) const { Erreurs::messageErreur("rempliTamponPentes non implemente pour melange utilise"); };
			virtual void recupereTamponPentes(double *tampon, int &compteur) { Erreurs::messageErreur("recupereTamponPentes non implemente pour melange utilise"); };

      //Accesseurs
      virtual double getDensite() const { Erreurs::messageErreur("getDensite non accessible pour le type de melange demande"); return 0.; };
      virtual double getPression() const { Erreurs::messageErreur("getPression non accessible pour le type de melange demande"); return 0.; };
      virtual double getU() const { Erreurs::messageErreur("getU non accessible pour le type de melange demande"); return 0.; };
      virtual double getV() const { Erreurs::messageErreur("getV non accessible pour le type de melange demande"); return 0.; };
      virtual double getW() const { Erreurs::messageErreur("getW non accessible pour le type de melange demande"); return 0.; };
      virtual Coord getVitesse() const { Erreurs::messageErreur("getVitesse non accessible pour le type de melange demande"); return 0; };
      virtual double getEnergie() const { Erreurs::messageErreur("getEnergie non accessible pour le type de melange demande"); return 0.; };
      virtual double getEnergieTotale() const { Erreurs::messageErreur("getEnergieTotale non accessible pour le type de melange demande"); return 0.; };
      virtual double getVitesseSonFigee() const { Erreurs::messageErreur("getVitesseSonFigee non accessible pour le type de melange demande"); return 0.; };
      virtual double getVitesseSonWood() const { Erreurs::messageErreur("getVitesseSonWood non accessible pour le type de melange demande"); return 0.; };

      virtual void setPression(const double &p) { Erreurs::messageErreur("setPression non implemente pour melange utilise"); };
      virtual void setVitesse(const double &u, const double &v, const double &w) { Erreurs::messageErreur("setVitesse non implemente pour melange utilise"); };
      virtual void setVitesse(const Coord &vit) { Erreurs::messageErreur("setVitesse non implemente pour melange utilise"); };
      virtual void setU(const double &u) { Erreurs::messageErreur("setU non implemente pour melange utilise"); };
      virtual void setV(const double &v) { Erreurs::messageErreur("setV non implemente pour melange utilise"); };
      virtual void setW(const double &w) { Erreurs::messageErreur("setW non implemente pour melange utilise"); };
      virtual void setEnergieTotale(double &energieTotale) { Erreurs::messageErreur("setEnergieTotale non accessible pour le type de melange demande"); };

      //Operateurs
      virtual void changeSigne() { Erreurs::messageErreur("changeSigne non implemente pour melange utilise"); };
      virtual void multiplieEtAjoute(const Melange &pentesMelangeTemp, const double &coeff) { Erreurs::messageErreur("multiplieEtAjoute non implemente pour melange utilise"); };
      virtual void diviser(const double &coeff) { Erreurs::messageErreur("diviser non implemente pour melange utilise"); };

    protected:
    private:

};

extern double *GlMel_Yk;

#endif // MELANGE_H
