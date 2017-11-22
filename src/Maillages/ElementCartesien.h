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

#ifndef ELEMENTCARTESIEN_H
#define ELEMENTCARTESIEN_H

#include "Element.h"

class ElementCartesien : public Element
{
public:
  ElementCartesien();
  virtual ~ElementCartesien();

  virtual void setVolume(const double &volume);
  virtual void setLCFL(const double &lCFL);
  virtual void setPos(const double &X, const double &Y, const double &Z);
  virtual void setPos(const Coord &pos);
  virtual void setPosX(const double &X);
  virtual void setPosY(const double &Y);
  virtual void setPosZ(const double &Z);

  //Pour methode AMR
  virtual void creerElementEnfant();
  virtual Element* getElementEnfant(const int &numeroEnfant);
  virtual void finaliseElementsEnfants();

protected:
  //Attributs pour methode AMR
  std::vector<ElementCartesien*> m_elementsEnfants;      /*!< Vecteur d'elements enfants */

private:
};

#endif // ELEMENTCARTESIEN_H