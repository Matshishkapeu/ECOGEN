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

#ifndef MELKAPILA_H
#define MELKAPILA_H

#include <vector>
#include "../Melange.h"

class MelKapila : public Melange
{
    public:
      MelKapila();
      MelKapila(tinyxml2::XMLElement *etat, std::string nomFichier);
      virtual ~MelKapila();

      virtual void ecritMelange(std::ofstream &fluxFichier) const;
      virtual void alloueEtCopieMelange(Melange **melange);
      virtual void alloueYk(const int &nombrePhases);
      virtual void copieMelange(Melange &melange);
      virtual double calculDensite(const double *alphak, const double *rhok, const int &nombrePhases);
      virtual double calculPression(const double *alphak, const double *pk, const int &nombrePhases);
      virtual double calculEnergieInterne(const double *Yk, const double *ek, const int &nombrePhases);
      virtual double calculVitesseSonFigee(const double *Yk, const double *ck, const int &nombrePhases);
      
      virtual double calculTsat(const Eos *eosLiq, const Eos *eosVap, const double &pression, double *dTsat=0);

      virtual void calculGrandeursMelange(Phase **vecPhase, const int &nombrePhases);
      virtual void energieInterneVersEnergieTotale(std::vector<GrandeursPhysAdd*> &vecGPA);
      virtual void energieTotaleVersEnergieInterne(std::vector<GrandeursPhysAdd*> &vecGPA);

      virtual void calculEnergieCapillaire();
      virtual void calculEnergieTotale();

      virtual void projection(const Coord &normale, const Coord &tangente, const Coord &binormale);
      virtual void projectionRepereAbsolu(const Coord &normale, const Coord &tangente, const Coord &binormale);

      //Specifique a l ecriture des donnees
      virtual int getNombreScalaires() const { return 2; };
      virtual int getNombreVecteurs() const { return 1; };
      virtual double renvoieScalaire(const int &numVar) const;
      virtual Coord* renvoieVecteur(const int &numVar);
      virtual std::string renvoieNomScalaire(const int &numVar) const;
      virtual std::string renvoieNomVecteur(const int &numVar) const;

      //Specifique au parallele
      virtual int nombreVariablesATransmettre() const;
      virtual void rempliTampon(double *tampon, int &compteur) const;
      virtual void recupereTampon(double *tampon, int &compteur);

      //Specifique a l ordre 2
      virtual void calculPentesMelange(const Melange &pGauche, const Melange &pDroite, const double &distance);
      virtual void miseAZero();
      virtual void extrapole(const Melange &pente, const double &distance);
      virtual void limitePentes(const Melange &penteGauche, const Melange &penteDroite, Limiteur &limiteurGlobal);

			//Specifique a l'ordre 2 parallele
			virtual int nombrePentesATransmettre() const;
			virtual void rempliTamponPentes(double *tampon, int &compteur) const;
			virtual void recupereTamponPentes(double *tampon, int &compteur);

      //Accesseurs
      virtual double getDensite() const;
      virtual double getPression() const;
      virtual double getU() const;
      virtual double getV() const;
      virtual double getW() const;
      virtual Coord getVitesse() const;
      virtual double getEnergie() const;
      virtual double getEnergieTotale() const;
      virtual double getVitesseSonFigee() const;
      virtual double getVitesseSonWood() const;

      virtual void setPression(const double &p);
      virtual void setVitesse(const double &u, const double &v, const double &w);
      virtual void setVitesse(const Coord &vit);
      virtual void setU(const double &u);
      virtual void setV(const double &v);
      virtual void setW(const double &w);
      virtual void setEnergieTotale(double &energieTotale);

      //Operateurs
      virtual void changeSigne();
      virtual void multiplieEtAjoute(const Melange &pentesMelangeTemp, const double &coeff);
      virtual void diviser(const double &coeff);

    protected:
    private:
      double m_densite;
      double m_pression;
      Coord m_vitesse;
      double m_energie;
      double m_energieTotale;
      double m_vitesseSonFigee;
      double m_vitesseSonWood;
};

#endif // MELKAPILA_H
