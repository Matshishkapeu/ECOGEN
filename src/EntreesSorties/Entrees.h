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

#ifndef ENTREES_H
#define ENTREES_H

#include <sstream>
#include <cassert>

#include "../libTierces/tinyxml2.h"
#include "../Eos/EnteteEquationEtat.h"
#include "../Geometries/EnteteDomaineGeometrique.h"
#include "../Sources/EnteteSource.h"
#include "../Ordre2/EnteteLimiteur.h"
#include "../CondLims/EnteteCondLim.h"
#include "../PhysiqueAdditionnelle/EnteteGPA.h"
#include "../Erreurs.h"
#include "../Outils.h"

class Entrees;

#include "../Run.h"
#include "Sorties.h"

class Entrees
{
  public:
    Entrees(Run *run);
    virtual ~Entrees();

    void lectureEntreesXML(std::vector<DomaineGeometrique*> &domaines, std::vector<CondLim*> &condLim);

    void entreeMain(std::string casTest);
    void entreeMaillage(std::string casTest);
    void entreeModele(std::string casTest);
    Eos* entreeEOS(std::string EOS, int &numeroEOS);
    void entreeConditionsInitiales(std::string casTest, std::vector<DomaineGeometrique*> &domaines, std::vector<CondLim*> &condLim);

	//Accesseur
	std::string getMain() const { return m_nomMain; };
	std::string getMaillage() const { return m_nomMaillage; };
	std::string getCI() const { return m_nomCI; };
	std::string getModele() const { return m_nomModele; };
  Run *getRun() const { return m_run; }

private:
  Run *m_run;    //pointeur vers run

	int m_vMain;		//!< Numero de version fichier entree main
	int m_vMaillage;    //!< Numero de version fichier entree maillage
	int m_vCI;          //!< Numero de version fichier entree conditionsInitiales
	int m_vModele;      //!< Numero de version fichier entree modele

	std::string m_nomMain;		//!< Nom du fichier main
	std::string m_nomMaillage;   //!< Nom du fichier maillage
	std::string m_nomCI;         //!< Nom du fichier conditionsInitiales
	std::string m_nomModele;     //!< Nom du fichier modele
};

#endif // ENTREES_H