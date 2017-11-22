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

#ifndef CELLULEO2GHOST_H
#define CELLULEO2GHOST_H

#include "CelluleO2.h"

class CelluleO2Ghost : public CelluleO2
{
public:
	CelluleO2Ghost();
	CelluleO2Ghost(int lvl); //Pour AMR
	virtual ~CelluleO2Ghost();

	virtual void alloue(const int &nombrePhases, const int &nombreTransports, const std::vector<PhysAdd*> &physAdd, Modele *modele);
	virtual void calculPentesLocal(const int &nombrePhases, const int &nombreTransports, BordDeMaille &bordRef, Limiteur &limiteurGlobal, Limiteur &limiteurInterface);
	virtual void creerCelluleEnfant(const int &num, const int &lvl);
	virtual void recupereTamponPentes(double *tampon, int &compteur);
	virtual void recupereTamponPentesAMR(double *tampon, int &compteur, const int &lvl);
	virtual bool estCelluleO2Ghost() const { return true; };

protected:
	Phase **m_vecPhasesPentesGhost;         /*!< vecteur des pentes des phases */
	Melange *m_melangePentesGhost;          /*!< vecteur des pentes de melange */
	double *m_vecTransportsPentesGhost;	    /*!< vecteur des pentes des transports */
	
private:
};

#endif // CELLULEO2GHOST_H
