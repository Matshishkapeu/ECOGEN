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

#include "OGDroite.h"
#include <iostream>

//***********************************************************************

OGDroite::OGDroite(){}

//***********************************************************************

OGDroite::OGDroite(const Coord &point, const Coord &vecDir) :
  ObjetGeometrique(DROITE), m_point(point), m_vecDir(vecDir)
{
  if (vecDir.norme() < 1e-6) {
    Erreurs::messageErreur("OGDroite::OGDroite impossible de creer droite, vecteur directeur nul");
  }
  else {
    m_vecDir.normalise();
  }
}

//***********************************************************************

OGDroite::~OGDroite(){}

//***********************************************************************

double OGDroite::distancePoint(const Coord &point) const
{
  Coord vec; vec.creeVecteur(m_point, point);
  if (fabs(vec.scalaire(m_vecDir)) < 1e-6) { return 0.; }
  else { return  (vec.vectoriel(m_vecDir)).norme(); }
}

//***********************************************************************

Coord OGDroite::projettePoint(const Coord &point) const
{
  Coord projete;
  Coord vec; vec.creeVecteur(m_point, point);
  projete.setXYZ(m_vecDir.scalaire(vec),0.,0.);
  return projete;
}

//***********************************************************************