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

#ifndef EOS_H
#define EOS_H
#include <string>
#include <vector>
#include <cassert>
#include "../Erreurs.h"
#include "../libTierces/tinyxml2.h"

class Eos
{
    public:
      Eos();
      Eos(int &numero);
      Eos(std::string nom, int &numero);
      virtual ~Eos();
      //Methode constantes
      void visualise() const;
      std::string retourneNom() const;
      int getNumero() const;

      //Lecture parametres physiques (viscosite, conductivite, etc.)
      void lectureParamPhysiques(tinyxml2::XMLNode *element, std::string nomFichier = "Fichier Inconnu");

      //Methodes generales pour toutes les EOS
      double calculEnthalpieTotale(const double &densite, const double &pression, const double &vitesse) const;
  
      //Methodes virtuelles pour les classes filles
      virtual void attributParametresEos(std::string nom, std::vector<double> parametresEos) = 0;

      virtual double calculTemperature(const double &densite, const double &pression) const=0; //Methodes virtuelles pures
      virtual double calculEnergie(const double &densite, const double &pression) const=0;
      virtual double calculPression(const double &densite, const double &energie) const=0;
      virtual double calculVitesseSon(const double &densite, const double &pression) const=0;
  
      virtual double calculPressionIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &densiteFinale) const=0;
      virtual double calculPressionHugoniot(const double &pressionInitiale, const double &densiteInitiale, const double &densiteFinale) const=0;
      virtual double calculDensiteIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *drhodp=0) const=0;
      virtual double calculDensiteHugoniot(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *drhodp = 0) const = 0;
      virtual double calculEnthalpieIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *dhdp = 0) const = 0;

      virtual double calculDensiteSaturation(const double &pression, const double &Tsat, const double &dTsatdP, double *drhodp = 0) const { Erreurs::messageErreur("calculDensiteSaturation non prevu pour EOS : " + m_nom); return 0; };
      virtual double calculRhoEnergieSaturation(const double &pression, const double &rho, const double &drhodp, double *drhoedp = 0) const { Erreurs::messageErreur("calculRhoEnergieSaturation non prevu pour EOS : " + m_nom); return 0; };

      virtual void renvoiSpecialEosMelange(double &gamPinfSurGamMoinsUn, double &eRef, double &unSurGamMoinsUn) const = 0;

      virtual double vfpfh(const double &pression, const double &enthalpie) const { Erreurs::messageErreur("vfpfh non prevu pour EOS : " + m_nom); return 0; };

      //Derivees partielles
      virtual double dvdpch(const double &pression, const double &enthalpie) const { Erreurs::messageErreur("dvdpch non prevu pour EOS : " + m_nom); return 0; };
      virtual double dvdhcp(const double &pression, const double &enthalpie) const { Erreurs::messageErreur("dvdhcp non prevu pour EOS : " + m_nom); return 0; };

      //Verifications
      virtual void verifiePression(const double &pression) const { Erreurs::messageErreur("verifiePression non prevu pour Eos : " + m_nom); };
      virtual void verifieEtCorrigePression(double &pression) const { Erreurs::messageErreur("verifieEtCorrigePression non prevu pour Eos : "+ m_nom); };

      //Mod
      virtual void renvoiInfo(double *&data) const = 0;

      //Accesseurs
      virtual double getGamma() const { return 0; };
      virtual double getPInf() const { return 0; };
      virtual double getCv() const { return 0; };
      virtual double getERef() const { return 0; };
      virtual double getSRef() const { return 0; };
      double getMu() const;
      double getLambda() const;

    protected:
      int m_numero;
      std::string m_nom;

      double m_mu;      //!< viscosite dynamique
      double m_lambda;  //!< conductivite thermique
};

#endif // EOS_H
