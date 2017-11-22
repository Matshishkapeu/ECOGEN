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

#ifndef OUTILS_H
#define OUTILS_H

#include <string>
#include "Eos/Eos.h"
#include <cmath>

typedef enum Axe{ X, Y, Z } Axe;
typedef enum TypeM { REC, UNS, AMR } TypeM;
//Enumeration pour le type de maillage (REC : rectilinear, UNS : unstructured, AMR : adaptative mesh refinement)
typedef enum TypeDonnee { FLOAT, DOUBLE, INT, CHAR } TypeDonnee;
typedef enum TypeOG { DROITE, PLAN } TypeOG;

class Outils
{
  public:
    Outils();
    Outils(const int &nombrePhases);
    virtual ~Outils();

    static void majuscule(std::string &chaine);
    static double pi();
    void alloueOutil(const int &nombrePhases = 1);

  public:

    double m_nombrePhases;

    double* ak;
    double* rhok;
    double* pk;
    double* akS;
    double* rhokS;
    Eos** eos;

};

extern Outils *BO;

#endif // OUTILS_H