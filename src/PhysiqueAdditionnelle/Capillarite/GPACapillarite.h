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

#ifndef GPACAPILLARITE_H
#define GPACAPILLARITE_H

#include "../GrandeursPhysAdd.h"

class GPACapillarite : public GrandeursPhysAdd
{
    public:
    GPACapillarite();
    GPACapillarite(PhysAdd* physadd);
    virtual ~GPACapillarite();

    virtual void calculGrandeurs(Cellule* cell);

    //Accesseurs
    virtual void setGrad(const Coord &grad, int num = -1);
    virtual Coord getGrad(int num = -1) const;

    protected:
    Coord m_gradC;                      /*!< gradient de la fonction de transport (vecteur w)*/

    private:
};

#endif // GPACAPILLARITE_H
