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

#ifndef ELEMENT_H
#define ELEMENT_H

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../Maths/Coord.h"
#include "../Maths/ObjetGeometrique.h"
#include "../Erreurs.h"
#include "../Outils.h"

class Element;

#include "Face.h"

class Element
{
public:
  Element();
  virtual ~Element();

  //Accesseurs
  void setCelluleAssociee(const int &numCellule);
  Coord getPosition() const;
  double getLCFL() const;
  double getVolume() const;
  int getNumCelluleAssociee() const;

  virtual int getIndice() const { Erreurs::messageErreur("getIndice inexistant pour Element demande"); return 0; };
  virtual void setVolume(const double &volume){ Erreurs::messageErreur("setVolume impossible pour element utilise"); };
  virtual void setLCFL(const double &lCFL){ Erreurs::messageErreur("setlCFL impossible pour element utilise"); };
  virtual void setPos(const double &X, const double &Y, const double &Z){ Erreurs::messageErreur("setPos impossible pour element utilise"); };
  virtual void setPos(const Coord &pos) { Erreurs::messageErreur("setPos impossible pour element utilise"); };
  virtual void setPosX(const double &X) { Erreurs::messageErreur("setPosX impossible pour element utilise"); };
  virtual void setPosY(const double &Y) { Erreurs::messageErreur("setPosY impossible pour element utilise"); };
  virtual void setPosZ(const double &Z) { Erreurs::messageErreur("setPosZ impossible pour element utilise"); };
  
  void ecritPos(std::ofstream &fluxFichier, Axe axe);

  virtual void afficheInfos() const{ Erreurs::messageErreur("AfficheInfos non prevue pour element utilise"); };
  
  Coord vecteur(const Element *e); /*!< Cree un vecteur a partir des centres d elements */
  Coord vecteur(const Face *f);    /*!< Cree un vecteur entre centre element et centre d une face */

  double distance(const Element *e);  /*!< Calcul de la distance entre centre et centre d un autre element */
  double distanceX(const Element *e); /*!< Calcul de la distance selon x entre centre et centre d un autre element */
  double distanceY(const Element *e); /*!< Calcul de la distance selon y entre centre et centre d un autre element */
  double distanceZ(const Element *e); /*!< Calcul de la distance selon z entre centre et centre d un autre element */
  double distance(const Face *f);     /*!< Calcul de la distance entre centre et centre d une face */
  double distanceX(const Face *f);    /*!< Calcul de la distance selon x entre centre et centre d une face */
  double distanceY(const Face *f);    /*!< Calcul de la distance selon y entre centre et centre d une face */
  double distanceZ(const Face *f);    /*!< Calcul de la distance selon z entre centre et centre d une face */

  bool traverseObjet(const ObjetGeometrique &objet) const;

  //Pour methode AMR
  virtual void creerElementEnfant() { Erreurs::messageErreur("creerElementsEnfants impossible pour element utilise"); };
  virtual Element* getElementEnfant(const int &numeroEnfant) { Erreurs::messageErreur("getElementEnfant impossible pour element utilise"); return 0; };
  virtual void finaliseElementsEnfants() { Erreurs::messageErreur("finaliseElementsEnfants impossible pour element utilise"); };

protected:

  Coord m_position;       /*!< Position du centre de l'element*/
  double m_volume;        /*!< Volume pour les elements 3D, Aire pour les elements 2D, Longueur pour les elements 1D*/
  double m_lCFL;          /*!< Longueur utile pour le calcul du pas de temps*/
  int m_numCelluleAssociee;
};

#endif // ELEMENT_H