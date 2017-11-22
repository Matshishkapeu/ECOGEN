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

#ifndef SORTIES_H
#define SORTIES_H

//Macro pour les interactions systeme (creation/destruction repertoires)
#ifdef WIN32
  #include <direct.h>
#else
  #include <sys/types.h>
  #include <sys/stat.h>
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "../libTierces/tinyxml2.h"
#include "../Erreurs.h"
#include "../Maillages/EnteteMaillage.h"
#include "../Cellule.h"
#include "IO.h"

class Sorties;

#include "Entrees.h"

class Sorties
{
  public:
    Sorties();
    Sorties(std::string casTest, std::string nomRun, tinyxml2::XMLElement *element, std::string nomFichier, Entrees *entree);
    virtual ~Sorties();

    void prepareSorties(const Cellule &cell);
    virtual void prepareSortiesInfos();
    void ecritSolution(Maillage *maillage, std::vector<Cellule *> *cellulesLvl);
    virtual void ecritInfos();

    virtual void prepareSortieSpecifique() { try { throw ErreurECOGEN("prepareSortieSpecifique non prevu pour sortie consideree"); } catch (ErreurECOGEN &) { throw; } };
    virtual void ecritSolutionSpecifique(Maillage *maillage, std::vector<Cellule *> *cellulesLvl) { try { throw ErreurECOGEN("ecritSolutionSpecifique non prevu pour sortie consideree"); } catch (ErreurECOGEN &) { throw; } };

    //Accesseur
    int getNumSortie() const { return m_numFichier; };

  protected:

    //Donnees generales
    void afficheInfoEcriture() const;
    void sauvegardeInfos() const;
    void sauvegardeInfosMailles() const;
    std::string creationNomFichier(const char* nom, int lvl = -1, int proc = -1, int numFichier = -1) const;

    void ecritJeuDonnees(std::vector<double> jeuDonnees, std::ofstream &fluxFichier, TypeDonnee typeDonnee);
    
	  Entrees *m_entree;												   //!<Pointeur vers entree
    Run *m_run;                                  //!<pointeur vers run

    //attribut nom fichiers/dossiers
    std::string m_casTest;                                             //!<Nom du cas test (defini dans "main.xml")
    std::string m_infosCalcul;                                         //!<Nom fichier pour sauvegarder les infos utiles du calcul
    std::string m_infoMailles;                                         //!<Nom fichiers pour stocker les infos de maillage
    std::string m_nomFichierResultats;                                 //!<Nom du fichier de sortie resultat
    std::string m_nomFichierCollection;                                //!<Nom de la collection regroupant les fichiers resultats
    std::string m_dossierSortie;                                       //!<Dossier pour enregistrement des resultats
    std::string m_dossierSauvegardesEntrees;                           //!<Dossier pour copier les fichiers entrees
    std::string m_dossierSauvegardesInfosMailles;                      //!<Dossier pour stocker les infos de maillage
    std::string m_dossierCoupes;
    std::string m_fichierCollection;                                   //!<Chemin du fichier collection regroupant les fichiers resultats
     
    //attribut parametres d ecriture
    bool m_ecritBinaire;                                               //!<Choix ecriture binaire/ASCII
    bool m_donneesSeparees;                                            //!<Choix ecriture donnees dans des fichiers separes

    int m_numFichier; 
    std::string m_endianMode;
    
    //Utile pour ecriture des donnees de cellules
    Cellule m_celluleRef;                                              //!<cellule de reference pour recupérer les noms des variables
};

#endif //SORTIES_H
