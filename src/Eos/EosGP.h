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

#ifndef EOSGP_H
#define EOSGP_H
#include "Eos.h"

class EosGP : public Eos
{
    public:
        EosGP();
        EosGP(std::vector<std::string> &nomParametresEos, int& numero);
        EosGP(std::string nom, double gamma, double cv, double eRef, int &numero);
        virtual ~EosGP();

        virtual void attributParametresEos(std::string nom, std::vector<double> parametresEos);

        //Methodes constantes (virtuelles car heritee de la classe Eos)
        virtual double calculTemperature(const double &densite,const double &pression) const;
        virtual double calculEnergie(const double &densite,const double &pression) const;
        virtual double calculPression(const double &densite,const double &energie) const;
        virtual double calculVitesseSon(const double &densite,const double &pression) const;

        virtual double calculPressionIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &densiteFinale) const;
        virtual double calculPressionHugoniot(const double &pressionInitiale, const double &densiteInitiale, const double &densiteFinale) const;
        virtual double calculDensiteIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *drhodp=0) const;
        virtual double calculDensiteHugoniot(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *drhodp = 0) const;
        virtual double calculEnthalpieIsentrope(const double &pressionInitiale, const double &densiteInitiale, const double &pressionFinale, double *dhdp = 0) const;

        virtual double calculDensiteSaturation(const double &pression, const double &Tsat, const double &dTsatdP, double *drhodp = 0) const;
        virtual double calculRhoEnergieSaturation(const double &pression, const double &rho, const double &drhodp, double *drhoedp = 0) const;

        virtual void renvoiSpecialEosMelange(double &gamPinfSurGamMoinsUn, double &eRef, double &unSurGamMoinsUn) const;

        virtual double vfpfh(const double &pression, const double &enthalpie) const;

        //derivees partielles
        virtual double dvdpch(const double &pression, const double &enthalpie) const;
        virtual double dvdhcp(const double &pression, const double &enthalpie) const;

        //Verifications
        virtual void verifiePression(const double &pression) const;
        virtual void verifieEtCorrigePression(double &pression) const;

		    //Mod
		    virtual void renvoiInfo(double *&data) const;

        //Accesseurs
        virtual double getGamma() const;
        virtual double getCv() const;
        virtual double getERef() const;
        virtual double getSRef() const;

    protected:
    private:
        double m_gamma;
        double m_cv;
        double m_eRef;
        double m_sRef;
};

#endif // EOSGP_H
