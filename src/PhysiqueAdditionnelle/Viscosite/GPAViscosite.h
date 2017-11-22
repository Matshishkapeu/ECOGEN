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

#ifndef GPAVISCOSITE_H
#define GPAVISCOSITE_H

#include "../GrandeursPhysAdd.h"

class GPAViscosite : public GrandeursPhysAdd
{
    public:
    GPAViscosite();
    GPAViscosite(PhysAdd* physAdd);
    virtual ~GPAViscosite();

    virtual void calculGrandeurs(Cellule* cell);

    //Accesseurs
    virtual void setGrad(const Coord &grad, int num = -1); //1:U, 2:V, 3:W
    virtual Coord getGrad(int num = -1) const;             //1:U, 2:V, 3:W

    protected:
    Coord m_gradU;      /*!< 1. vecteur du gradient de la vitesse selon x de la cellule*/
    Coord m_gradV;      /*!< 2. vecteur du gradient de la vitesse selon y de la cellule*/
    Coord m_gradW;      /*!< 3. vecteur du gradient de la vitesse selon z de la cellule*/

    private:
};

#endif // GPAVISCOSITE_H
