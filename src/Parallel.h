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

#ifndef PARALLEL_H
#define PARALLEL_H

#include <mpi.h>
#include "Modeles/Phase.h"
#include "Cellule.h"

class Parallel
{
public:
  Parallel();
  virtual ~Parallel();

  void initialisation(int &argc, char *argv[]);
  void setVoisin(const int voisin); //attribut voisin
  void setElementsAEnvoyer(const int voisin, int* numeroElement, const int &nombreElements);
  void setElementsARecevoir(const int voisin, int* numeroElement, const int &nombreElements);
	void initialiseCommunicationsPersistantes(const int &nombreVariablesPrimitives, const int &nombreVariablesPentes, const int &nombreVariablesTransports, const int &dim);
	void calculDt(double &dt);
	void finalise(const int &lvlMax);
	void arreterCode();
	void verifieEtatCPUs();
  
  //Methodes pour toutes les variables primitives
  void initialiseCommunicationsPersistantesPrimitives();
  void finaliseCommunicationsPersistantesPrimitives(const int &lvlMax);
  void communicationsPrimitives(Cellule **cellules, Eos **eos, Prim type = vecPhases);

	//Methodes pour toutes les pentes
	void initialiseCommunicationsPersistantesPentes();
	void finaliseCommunicationsPersistantesPentes(const int &lvlMax);
	void communicationsPentes(Cellule **cellules);

  //Methodes pour une variable scalaire
  void initialiseCommunicationsPersistantesScalaire();
  void finaliseCommunicationsPersistantesScalaire(const int &lvlMax);
  void communicationsScalaire(Cellule **cellules, std::string nomScalaire);
  
  //Methodes pour une variable vectorielle
  void initialiseCommunicationsPersistantesVecteur(const int &dim);
  void finaliseCommunicationsPersistantesVecteur(const int &lvlMax);
  void communicationsVecteur(Cellule **cellules, std::string nomVecteur, const int &dim, int num=0, int indice=-1);

  //Methodes pour toutes les variables primitives
  void initialiseCommunicationsPersistantesTransports();
  void finaliseCommunicationsPersistantesTransports(const int &lvlMax);
  void communicationsTransports(Cellule **cellules);

	//Methodes pour les variables AMR
	void initialiseCommunicationsPersistantesAMR(const int &nombreVariablesPrimitives, const int &nombreVariablesPentes, const int &nombreVariablesTransports, const int &dim, const int &lvlMax);
	void initialiseCommunicationsPersistantesNiveauxAMR(const int &lvlMax);
	void miseAJourCommunicationsPersistantesLvl(int lvl, const int &dim);
	void finaliseAMR(const int &lvlMax);

	void initialiseCommunicationsPersistantesXi();
	void finaliseCommunicationsPersistantesXi(const int &lvlMax);
	void communicationsXi(Cellule **cellules, const int &lvl);

	void initialiseCommunicationsPersistantesSplit();
	void finaliseCommunicationsPersistantesSplit(const int &lvlMax);
	void communicationsSplit(Cellule **cellules, const int &lvl);

	void initialiseCommunicationsPersistantesNombreCellulesGhost();
	void finaliseCommunicationsPersistantesNombreCellulesGhost();
	void communicationsNombreCellulesGhost(Cellule **cellules, const int &lvl);

	void communicationsPrimitivesAMR(Cellule **cellules, Eos **eos, const int &lvl, Prim type = vecPhases);
	void communicationsPentesAMR(Cellule **cellules, const int &lvl);
	void communicationsScalaireAMR(Cellule **cellules, std::string nomScalaire, const int &lvl);
	void communicationsVecteurAMR(Cellule **cellules, std::string nomVecteur, const int &dim, const int &lvl, int num = 0, int indice = -1);
  void communicationsTransportsAMR(Cellule **cellules, const int &lvl);

private:
    
  int m_etatCPU;
  bool *m_estVoisin;
	int ** m_elementsAEnvoyer;
	int ** m_elementsARecevoir;
	int * m_nombreElementsAEnvoyerAVoisin;
	int * m_nombreElementsARecevoirDeVoisin;
	int m_nombreVariablesPrimitives;          /*Nombre de variables de primitives a transmettre (phases + melange + transports)*/
	int m_nombreVariablesPentes;              /*Nombre de variables de pentes a transmettre (phases + melange + transports)*/
  int m_nombreVariablesTransports;          /*Nombre de variables de transports a transmettre*/

	std::vector<double **> m_tamponRec;
	std::vector<double **> m_tamponEnv;
	std::vector<double **> m_tamponRecPentes;
	std::vector<double **> m_tamponEnvPentes;
	std::vector<double **> m_tamponRecScalaire;
	std::vector<double **> m_tamponEnvScalaire;
	std::vector<double **> m_tamponRecVecteur;
	std::vector<double **> m_tamponEnvVecteur;
  std::vector<double **> m_tamponRecTransports;
  std::vector<double **> m_tamponEnvTransports;
	std::vector<double **> m_tamponRecXi;
	std::vector<double **> m_tamponEnvXi;
	std::vector<bool **> m_tamponRecSplit;
	std::vector<bool **> m_tamponEnvSplit;
	int * m_tamponNombreElementsAEnvoyerAVoisin;
	int * m_tamponNombreElementsARecevoirDeVoisin;
  
	std::vector<MPI_Request **> m_reqEnvois;
	std::vector<MPI_Request **> m_reqReceptions;
	std::vector<MPI_Request **> m_reqEnvoisPentes;
	std::vector<MPI_Request **> m_reqReceptionsPentes;
	std::vector<MPI_Request **> m_reqEnvoisScalaire;
	std::vector<MPI_Request **> m_reqReceptionsScalaire;
	std::vector<MPI_Request **> m_reqEnvoisVecteur;
	std::vector<MPI_Request **> m_reqReceptionsVecteur;
  std::vector<MPI_Request **> m_reqEnvoisTransports;
  std::vector<MPI_Request **> m_reqReceptionsTransports;
	std::vector<MPI_Request **> m_reqEnvoisXi;
	std::vector<MPI_Request **> m_reqReceptionsXi;
	std::vector<MPI_Request **> m_reqEnvoisSplit;
	std::vector<MPI_Request **> m_reqReceptionsSplit;
	MPI_Request ** m_reqNombreElementsAEnvoyerAVoisin;
	MPI_Request ** m_reqNombreElementsARecevoirDeVoisin;

};

extern Parallel Calcul_Parallele;
extern int rang;
extern int Ncpu;

#endif // PARALLEL_H