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

#include "PhysAdd.h"

using namespace std;

//***********************************************************************

PhysAdd::PhysAdd(){}

//***********************************************************************

PhysAdd::~PhysAdd(){}

//***********************************************************************

void PhysAdd::calculFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases)
{
  this->resolFluxPhysAdd(bord, nombrePhases);

  if (bord->getCellGauche()->getLvl() == bord->getCellDroite()->getLvl()) {     //CoefAMR = 1 pour les deux
    this->ajoutFluxPhysAdd(bord, nombrePhases, 1.);                             //Ajout du flux sur maille droite
    this->retireFluxPhysAdd(bord, nombrePhases, 1.);                            //Retrait du flux sur maille gauche
  }
  else if (bord->getCellGauche()->getLvl() > bord->getCellDroite()->getLvl()) { //CoefAMR = 1 pour la gauche et 0.5 pour la droite
    this->ajoutFluxPhysAdd(bord, nombrePhases, 0.5);                            //Ajout du flux sur maille droite
    this->retireFluxPhysAdd(bord, nombrePhases, 1.);                            //Retrait du flux sur maille gauche
  }
  else {                                                                        //CoefAMR = 0.5 pour la gauche et 1 pour la droite
    this->ajoutFluxPhysAdd(bord, nombrePhases, 1.);                             //Ajout du flux sur maille droite
    this->retireFluxPhysAdd(bord, nombrePhases, 0.5);                           //Retrait du flux sur maille gauche
  }
}

//***********************************************************************

void PhysAdd::ajoutFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases, const double &coefAMR)
{
  double volume(bord->getCellDroite()->getElement()->getVolume());
  double surface(bord->getFace()->getSurface());
  double coefA(surface / volume); //pas de "pas de temps"
  coefA = coefA*coefAMR;
  bord->getCellDroite()->getCons()->ajoutFlux(coefA, nombrePhases);
}

//***********************************************************************

void PhysAdd::retireFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases, const double &coefAMR)
{
  double volume(bord->getCellGauche()->getElement()->getVolume());
  double surface(bord->getFace()->getSurface());
  double coefA(surface / volume); //pas de "pas de temps"
  coefA = coefA*coefAMR;
  bord->getCellGauche()->getCons()->retireFlux(coefA, nombrePhases);
}

//***********************************************************************

void PhysAdd::calculFluxPhysAddLimite(BordDeMaille *bord, const int &nombrePhases)
{
  this->resolFluxPhysAddLimite(bord, nombrePhases);
  this->retireFluxPhysAdd(bord, nombrePhases, 1.); //Retrait du flux sur maille gauche
}

//***********************************************************************

void PhysAdd::ajoutNonConsPhysAdd(Cellule *cell, const int &nombrePhases)
{
  this->ajoutNonCons(cell, nombrePhases);
}

//***********************************************************************
