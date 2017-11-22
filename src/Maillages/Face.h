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

#ifndef FACE_H
#define FACE_H

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include "../Maths/Coord.h"
#include "../Erreurs.h"

class Face;

#include "Element.h"

class Face
{
public:
  Face();
  virtual ~Face();

  //Accesseurs
  Coord getNormale() const;
  Coord getTangente() const;
  Coord getBinormale() const;
  double getSurface() const;
  Coord getPos() const;

  virtual void setSurface(const double &surface){ Erreurs::messageErreur("setSurface non prevue pour le type de face demande"); };
  virtual void initialiseAutres(const double &surface, const Coord &normale, const Coord &tangente, const Coord &binormale){ Erreurs::messageErreur("initialiseAutres non prevue pour le type de face demande"); }
  virtual void setPos(const double &X, const double &Y, const double &Z) { Erreurs::messageErreur("setPos impossible pour face utilisee"); };
  virtual void setNormale(const double &X, const double &Y, const double &Z) { Erreurs::messageErreur("setNormale impossible pour face utilisee"); };
  virtual void setTangente(const double &X, const double &Y, const double &Z) { Erreurs::messageErreur("setTangente impossible pour face utilisee"); };
  virtual void setBinormale(const double &X, const double &Y, const double &Z) { Erreurs::messageErreur("setBinormale impossible pour face utilisee"); };

  Coord vecteur(Element *e);   /*!< Cree vecteur entre centre face et centre d un element */
  double distance(Element *e); /*!< Calcul de la distance a un centre d element */

  virtual void afficheInfos() const{ Erreurs::messageErreur("AfficheInfos non prevue pour le type de face demande"); };

  //Pour methode AMR
  virtual Face* creerNouvelleFace() { Erreurs::messageErreur("creerNouvelleFace impossible pour face utilise"); return 0; };

protected:

  Coord m_position;     /*!< Position du centre de la face */
  double m_surface;     /*!< 1.0 pour element 0D, longueur pour element 1D, surface pour element 2D */
  Coord m_normale;
  Coord m_tangente;
  Coord m_binormale;

};

#endif // FACE_H
