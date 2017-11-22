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

#ifndef GRANDEURSPHYSADD_H
#define GRANDEURSPHYSADD_H

class GrandeursPhysAdd;

#include "PhysAdd.h"

class GrandeursPhysAdd
{
    public:
    GrandeursPhysAdd();
    GrandeursPhysAdd(PhysAdd* physadd);
    virtual ~GrandeursPhysAdd();

    virtual void calculGrandeurs(Cellule* cell) { Erreurs::messageErreur("calculsGrandeurs non prevue pour Grandeur Physique additionelle"); };
    double calculEnergiePhysAdd();

    //Accesseurs pour les differents gradients
    virtual void setGrad(const Coord &grad, int num=-1) { Erreurs::messageErreur("setGrad non prevue pour Grandeur Physique additionelle"); };
    virtual Coord getGrad(int num=-1) const { Erreurs::messageErreur("getGrad non prevue pour Grandeur Physique additionelle"); return 0; };

    virtual void setGradU(const Coord &grad) {};
    virtual void setGradV(const Coord &grad) {};
    virtual void setGradW(const Coord &grad) {};
    virtual void setGradTk(int &numPhase, const Coord &grad) {};
    virtual Coord getGradU() const { return 0; };
    virtual Coord getGradV() const { return 0; };
    virtual Coord getGradW() const { return 0; };
    virtual Coord getGradTk(int &numPhase) const { return 0; };

    PhysAdd* getPhysAdd() { return m_physAdd; };

    protected:
      PhysAdd* m_physAdd;

    private:
};

#endif // GRANDEURSPHYSADD_H
