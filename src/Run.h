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

#ifndef RUN_H
#define RUN_H

class Run;

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <sstream>
#include "Cellule.h"
#include "Modeles/EntetePhase.h"
#include "BordDeMaille.h"
#include "Parallel.h"
#include "Maillages/EnteteMaillage.h"
#include "CondLims/EnteteCondLim.h"
#include "Eos/EnteteEquationEtat.h"
#include "Modeles/EnteteModele.h"
#include "Geometries/EnteteDomaineGeometrique.h"
#include "Ordre2/EnteteLimiteur.h"
#include "PhysiqueAdditionnelle/EnteteGPA.h"
#include "PhysiqueAdditionnelle/EntetePhysAdd.h"

#include "EntreesSorties/Entrees.h"
#include "EntreesSorties/Sorties.h"
//#include "EntreesSorties/EnteteEntreesSorties.h"

//FP//TODO// Mettre La classe globale en classe amie de pas mal d autres classes pour simplifier les appels de methodes

class Run
{
  public:
    /** Default constructor */
    Run(std::string nomCasTest, const int &numero);
    //Run(std::string);
    /** Default destructor */
    virtual ~Run();

    void initialisation(int argc, char* argv[]);
    void repriseFichier(int &iteration, double &dt, double &tempsPhysique, clock_t &temps);

    /*! \brief Methode de resolution Hyperbolique + Relaxations + Termes Sources
     *
     *  La partie hyperbolique est resolue sous la forme :
     *  \f[
     *       \frac{U^{n+1}_i-U^{n}_i}{\Delta t} =  -\sum_{faces} \vec{F}^*_f \cdot \vec{n}_f
     *  \f]
     */
    void resolution();
    void procedureIntegration(double &dt, int lvl, double &dtMax, int &nbMaillesTotalAMR);
    void procedureAvancement(double &dt, int &lvl, double &dtMax) const;
    void resolHyperbolique(double &dt, int &lvl, double &dtMax) const;
    void resolHyperboliqueO2(double &dt, int &lvl, double &dtMax) const;
    void resolPhysiquesAdditionelles(double &dt, int &lvl) const;
    void resolTermesSources(double &dt, int &lvl) const;
    void resolRelaxations(int &lvl) const;
    void verifieErreurs() const;
    void finalise();
    static void arretApresErreur();

    //Acceseur
    int getNombrePhases() const;

  private:   

    int m_numTest;                            //!<Numero du cas test

    //Attribut entrees
    std::string m_casTest;                     //!<Contient le nom du cas test
    bool m_controleIterations;
    int m_nbIte, m_freq;
    float m_tempsFinal, m_freqTemps;
    double m_cfl;
    int m_nombrePhases;
    int m_nombreTransports;
    int m_nombreEos;                           //!<Nombre d equations d etat
    std::string m_ordre;                       //!<Ordre du calcul (ordre1, ordre2)
    int m_nombrePhysAdd;
    int m_nombreSources;
    int m_dimension;                           //!<dimension 1, 2 ou 3
    int m_evaporation;                         //!<indicateur evaporation 0 ou 1

    //Pour methode AMR
    int m_lvlMax;                              //!<Niveau maximal sur l arbre AMR (si m_lvlMax = 0, pas d AMR)
    int m_nbMaillesTotalAMR;                   //!<Nombre de mailles total maximum durant la simulation
    std::vector<Cellule *> *m_cellulesLvl;     //!<Tableau de vecteurs contenant les cellules de calcul, un vecteur par niveau.
    std::vector<BordDeMaille *> *m_bordsLvl;   //!<Tableau de vecteurs contenant les bords de calcul, un vecteur par niveau.

    //Attribut geometrie
    bool m_pretraitementParallele;             //!<variable sur le choix du pretraitement parallele (necessaire si la geometrie n a jamais ete cree)
    
    //Attribut Calcul
    Maillage *m_maillage;                      //!<Objet contenant toutes les proprietes geometrique des elements de maillage
    Modele *m_modele;                          //!<Objet contenant le modele physique
    Cellule **m_cellules;                      //!<Objet contenant les cellules de calcul (variables des phases, etc.)
    BordDeMaille **m_bords;                    //!<Tableau des faces entres cellules ou entre une cellule de calcul et une limite
    Eos **m_eos;                               //!<Tableau des equations d etats
    std::vector<PhysAdd*> m_physAdd;           //!<Vecteur des physiques additionnelles
    std::vector<Source*> m_sources;            //!<Vecteur de termes sources
    Limiteur *m_limiteurGlobal;                //!<Objet contenant le limiteur de pente pour ordre 2 en espace
    Limiteur *m_limiteurInterface;             //!<Objet contenant le limiteur de pente pour ordre 2 en espace sur l'interface (alpha, transports)
    std::vector<std::string> m_nomGTR;         //!<Vecteur des noms des grandeurs transportees 
    std::vector<std::string> m_nomGPA;         //!<Vecteur des noms des grandeurs physiques additionelles 
    std::vector<std::string> m_nomGPH;         //!<Vecteur des noms des grandeurs des phases
    double m_dt;
    double m_tempsPhys;
    clock_t m_tempsCalc;
    int m_iteration;

    int m_repriseFichier;                      //!<Numero du fichier si reprise de calcul

    //Attribut entrees / sorties
	  Entrees* m_entree;						               //!<Entree
    Sorties* m_sortie;                         //!<Sortie principale
    std::vector<Sorties *> m_coupes;           //!<Coupes

    friend class Entrees;
    friend class Sorties;
    friend class SortiesXML;
    friend class SortiesGNU;
    friend class Maillage;
    
    //Analyse temps - Les attributs sont stockes en millisecondes (a diviser par CLOCKS_PER_SEC)
    clock_t m_tempsInitial;

};

#endif // RUN_H
