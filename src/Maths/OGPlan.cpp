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

#include "OGPlan.h"

using namespace std;

//***********************************************************************

OGPlan::OGPlan(){}

//***********************************************************************

OGPlan::OGPlan(const Coord &point, const Coord &normale) :
  ObjetGeometrique(PLAN), m_point(point), m_normale(normale)
{
  if (normale.norme() < 1e-6) { 
    Erreurs::messageErreur("OGPlan::OGPlan impossible de creer plan, vecteur normal nul"); 
  }
  else {
    m_normale.normalise();
  }
  creeBase();
}


//***********************************************************************

OGPlan::~OGPlan(){}

//***********************************************************************

double OGPlan::distancePoint(const Coord &point) const
{
  Coord vec; vec.creeVecteur(m_point, point);
  return  fabs(vec.scalaire(m_normale));
}

//***********************************************************************

Coord OGPlan::projettePoint(const Coord &point) const
{
  Coord projete;
  Coord vec; vec.creeVecteur(m_point, point);
  projete.setXYZ(m_tangente.scalaire(vec), m_binormale.scalaire(vec), 0.);
  return projete;
}

//***********************************************************************

void OGPlan::creeBase() 
{
  Coord M, N;
  M = m_point;

  //Determination d'une tangente
  if (fabs(m_normale.getZ()) >= 1e-6) {
    N.setX(0); N.setY(0);
    N.setZ((M.getX()*m_normale.getX() + M.getY()*m_normale.getY()) / m_normale.getZ() + M.getZ());
  }
  else if (fabs(m_normale.getY()) >= 1e-6) {
    N.setX(0); N.setZ(0);
    N.setY((M.getX()*m_normale.getX() + M.getZ()*m_normale.getZ()) / m_normale.getY() + M.getY());
  }
  else if (fabs(m_normale.getX()) >= 1e-6) {
    N.setY(0); N.setZ(0);
    N.setX((M.getY()*m_normale.getY() + M.getZ()*m_normale.getZ()) / m_normale.getX() + M.getX());
  }
  else {
    Erreurs::messageErreur("OGPlan::creeBase impossible, vecteur normal problematique");
  }
  m_tangente.creeVecteur(M, N);
  m_tangente.normalise();
  m_binormale = Coord::produitVectoriel(m_normale, m_tangente);
}

//***********************************************************************