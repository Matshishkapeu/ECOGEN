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

#include "Parallel.h"
#include <iostream>
#include "Eos/Eos.h"

//Variables liees au calcul parallele
Parallel Calcul_Parallele;
int rang, Ncpu;

using namespace std;

//***********************************************************************

Parallel::Parallel() :
  m_etatCPU(1)
{}

//***********************************************************************

Parallel::~Parallel(){}

//***********************************************************************

void Parallel::initialisation(int &argc, char* argv[])
{
  //MPI_Init(&argc, &argv);
  //MPI_Comm_rank(MPI_COMM_WORLD, &rang);
  //MPI_Comm_size(MPI_COMM_WORLD, &Ncpu);
  if (Ncpu == 1) return; //Pas besoin de la suite dans le cas monoCPU

  m_estVoisin = new bool[Ncpu];
	m_elementsAEnvoyer = new int*[Ncpu];
  m_elementsARecevoir = new int*[Ncpu];
  m_nombreElementsAEnvoyerAVoisin = new int[Ncpu];
  m_nombreElementsARecevoirDeVoisin = new int[Ncpu]; // A priori peuvent etre differents si maillage tordu !

  m_tamponEnv.push_back(new double*[Ncpu]);
  m_tamponRec.push_back(new double*[Ncpu]);
	m_tamponEnvPentes.push_back(new double*[Ncpu]);
	m_tamponRecPentes.push_back(new double*[Ncpu]);
  m_tamponEnvScalaire.push_back(new double*[Ncpu]);
  m_tamponRecScalaire.push_back(new double*[Ncpu]);
  m_tamponEnvVecteur.push_back(new double*[Ncpu]);
  m_tamponRecVecteur.push_back(new double*[Ncpu]);
  m_tamponEnvTransports.push_back(new double*[Ncpu]);
  m_tamponRecTransports.push_back(new double*[Ncpu]);
	m_tamponEnvXi.push_back(new double*[Ncpu]);
	m_tamponRecXi.push_back(new double*[Ncpu]);
	m_tamponEnvSplit.push_back(new bool*[Ncpu]);
	m_tamponRecSplit.push_back(new bool*[Ncpu]);
	m_tamponNombreElementsAEnvoyerAVoisin = new int[Ncpu];
	m_tamponNombreElementsARecevoirDeVoisin = new int[Ncpu];

  m_reqEnvois.push_back(new MPI_Request*[Ncpu]);
  m_reqReceptions.push_back(new MPI_Request*[Ncpu]);
	m_reqEnvoisPentes.push_back(new MPI_Request*[Ncpu]);
	m_reqReceptionsPentes.push_back(new MPI_Request*[Ncpu]);
  m_reqEnvoisScalaire.push_back(new MPI_Request*[Ncpu]);
  m_reqReceptionsScalaire.push_back(new MPI_Request*[Ncpu]);
  m_reqEnvoisVecteur.push_back(new MPI_Request*[Ncpu]);
  m_reqReceptionsVecteur.push_back(new MPI_Request*[Ncpu]);
  m_reqEnvoisTransports.push_back(new MPI_Request*[Ncpu]);
  m_reqReceptionsTransports.push_back(new MPI_Request*[Ncpu]);
	m_reqEnvoisXi.push_back(new MPI_Request*[Ncpu]);
	m_reqReceptionsXi.push_back(new MPI_Request*[Ncpu]);
	m_reqEnvoisSplit.push_back(new MPI_Request*[Ncpu]);
	m_reqReceptionsSplit.push_back(new MPI_Request*[Ncpu]);
	m_reqNombreElementsAEnvoyerAVoisin = new MPI_Request*[Ncpu];
	m_reqNombreElementsARecevoirDeVoisin = new MPI_Request*[Ncpu];

  for (int i = 0; i < Ncpu; i++) {
    m_estVoisin[i] = false;
    m_nombreElementsAEnvoyerAVoisin[i] = 0;
    m_nombreElementsARecevoirDeVoisin[i] = 0;
		m_elementsAEnvoyer[i] = NULL;
    m_elementsARecevoir[i] = NULL;
    m_tamponEnv[0][i] = NULL;
    m_tamponRec[0][i] = NULL;
    m_reqEnvois[0][i] = NULL;
    m_reqReceptions[0][i] = NULL;
		m_tamponEnvPentes[0][i] = NULL;
		m_tamponRecPentes[0][i] = NULL;
		m_reqEnvoisPentes[0][i] = NULL;
		m_reqReceptionsPentes[0][i] = NULL;
    m_tamponEnvScalaire[0][i] = NULL;
    m_tamponRecScalaire[0][i] = NULL;
    m_reqEnvoisScalaire[0][i] = NULL;
    m_reqReceptionsScalaire[0][i] = NULL;
    m_tamponEnvVecteur[0][i] = NULL;
    m_tamponRecVecteur[0][i] = NULL;
    m_reqEnvoisVecteur[0][i] = NULL;
    m_reqReceptionsVecteur[0][i] = NULL;
    m_tamponEnvTransports[0][i] = NULL;
    m_tamponRecTransports[0][i] = NULL;
    m_reqEnvoisTransports[0][i] = NULL;
    m_reqReceptionsTransports[0][i] = NULL;
		m_tamponEnvXi[0][i] = NULL;
		m_tamponRecXi[0][i] = NULL;
		m_reqEnvoisXi[0][i] = NULL;
		m_reqReceptionsXi[0][i] = NULL;
		m_tamponEnvSplit[0][i] = NULL;
		m_tamponRecSplit[0][i] = NULL;
		m_reqEnvoisSplit[0][i] = NULL;
		m_reqReceptionsSplit[0][i] = NULL;
		m_tamponNombreElementsAEnvoyerAVoisin[i] = 0;
		m_tamponNombreElementsARecevoirDeVoisin[i] = 0;
		m_reqNombreElementsAEnvoyerAVoisin[i] = NULL;
		m_reqNombreElementsARecevoirDeVoisin[i] = NULL;
  } 
}

//***********************************************************************

void Parallel::setVoisin(const int voisin)
{ 
  m_estVoisin[voisin] = true;
}

//***********************************************************************

void Parallel::setElementsAEnvoyer(const int voisin, int* numeroElement, const int &nombreElements)
{
  m_nombreElementsAEnvoyerAVoisin[voisin] = nombreElements;

  //On dimensionne le tableau où seront stockes les numeros des mailles a envoyer au voisin "voisin"
	m_elementsAEnvoyer[voisin] = new int[nombreElements];
  //Puis on stocke les numeros a partir du tableau fourni
  for (int i = 0; i < nombreElements; i++) {
		m_elementsAEnvoyer[voisin][i] = numeroElement[i];
  }
}

//***********************************************************************

void Parallel::setElementsARecevoir(const int voisin, int* numeroElement, const int &nombreElements)
{
  m_nombreElementsARecevoirDeVoisin[voisin] = nombreElements;
  
  //On dimensionne le tableau où seront stockes les numeros des mailles où seront reçues les datas emises par le voisin "voisin"
  m_elementsARecevoir[voisin] = new int[nombreElements];
  //Puis on stocke les numeros a partir du tableau fourni
  for (int i = 0; i < nombreElements; i++) {
    m_elementsARecevoir[voisin][i] = numeroElement[i];
  }
}

//***********************************************************************

void Parallel::initialiseCommunicationsPersistantes(const int &nombreVariablesPrimitives, const int &nombreVariablesPentes, const int &nombreVariablesTransports, const int &dim)
{
	if (Ncpu > 1) {
		m_nombreVariablesPrimitives = nombreVariablesPrimitives;
    m_nombreVariablesPentes = nombreVariablesPentes;
    m_nombreVariablesTransports = nombreVariablesTransports;
		//Initialisation des communications des variables primitives du modele resolu
		Calcul_Parallele.initialiseCommunicationsPersistantesPrimitives();
		//Initialisation des communications des pentes pour l'ordre 2
		Calcul_Parallele.initialiseCommunicationsPersistantesPentes();
		//Initialisation des communications necessaires aux physiques additionelles (vecteurs de dim=3)
		Calcul_Parallele.initialiseCommunicationsPersistantesVecteur(dim);
    //Initialisation des communications des variables transportees
    Calcul_Parallele.initialiseCommunicationsPersistantesTransports();
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::calculDt(double &dt)
{
	double dt_temp = dt;
	MPI_Allreduce(&dt_temp, &dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::finalise(const int &lvlMax)
{
	if (Ncpu > 1) {
		this->finaliseCommunicationsPersistantesPrimitives(lvlMax);
		this->finaliseCommunicationsPersistantesPentes(lvlMax);
		this->finaliseCommunicationsPersistantesVecteur(lvlMax);
    this->finaliseCommunicationsPersistantesTransports(lvlMax);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::arreterCode()
{
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	exit(0);
}

//***********************************************************************

void Parallel::verifieEtatCPUs()
{
	//Regroupement des erreurs
	int nbErr_temp(0);
	int nbErr(erreurs.size());
	MPI_Allreduce(&nbErr, &nbErr_temp, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
	//Arret si erreur sur un CPU
	if (nbErr_temp) {
		Erreurs::arretCodeApresErreur(erreurs);
	}
}

//****************************************************************************
//************* Methodes pour toutes les variables primitives ****************
//****************************************************************************

void Parallel::initialiseCommunicationsPersistantesPrimitives()
{
  for (int voisin = 0; voisin < Ncpu; voisin++) {
    if (m_estVoisin[voisin]) { //Si CPU voisin
      //Determination du nombre de variables a communiquer
      int nombreEnvoi = m_nombreVariablesPrimitives*m_nombreElementsAEnvoyerAVoisin[voisin];
      int nombreRecoi = m_nombreVariablesPrimitives*m_nombreElementsARecevoirDeVoisin[voisin];

      //Nouvelle requête d'envois et son buffer associe
      m_reqEnvois[0][voisin] = new MPI_Request;
      m_tamponEnv[0][voisin] = new double[nombreEnvoi];
      MPI_Send_init(m_tamponEnv[0][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvois[0][voisin]);

      //Nouvelle requête de receptions et son buffer associe
      m_reqReceptions[0][voisin] = new MPI_Request;
      m_tamponRec[0][voisin] = new double[nombreRecoi];
      MPI_Recv_init(m_tamponRec[0][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptions[0][voisin]);
    }
  }
}

//***********************************************************************

void Parallel::finaliseCommunicationsPersistantesPrimitives(const int &lvlMax)
{
	for (int lvl = 0; lvl <= lvlMax; lvl++) {
		for (int voisin = 0; voisin < Ncpu; voisin++)	{
			if (m_estVoisin[voisin]) { //Si CPU voisin
				MPI_Request_free(m_reqEnvois[lvl][voisin]);
				MPI_Request_free(m_reqReceptions[lvl][voisin]);
        delete m_reqEnvois[lvl][voisin];
        delete[] m_tamponEnv[lvl][voisin];
        delete m_reqReceptions[lvl][voisin];
        delete[] m_tamponRec[lvl][voisin];
			}
		}
		delete[] m_reqEnvois[lvl];
		delete[] m_tamponEnv[lvl];
		delete[] m_reqReceptions[lvl];
		delete[] m_tamponRec[lvl];
	}
  m_reqEnvois.clear();
  m_tamponEnv.clear();
  m_reqReceptions.clear();
  m_tamponRec.clear();

}

//***********************************************************************

void Parallel::communicationsPrimitives(Cellule **cellules, Eos **eos, Prim type)
{
  int count(0);
  MPI_Status status;

  for (int voisin = 0; voisin < Ncpu; voisin++) {
    if (m_estVoisin[voisin]) {
      //Prepation des envois
      count = -1;
      for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
				cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponPrimitives(m_tamponEnv[0][voisin], count, type);
      }

      //Requête d'envoi
      MPI_Start(m_reqEnvois[0][voisin]);
      //Requête de reception
      MPI_Start(m_reqReceptions[0][voisin]);
      //Attente
      MPI_Wait(m_reqEnvois[0][voisin], &status);
      MPI_Wait(m_reqReceptions[0][voisin], &status);

      //Receptions
      count = -1;
      for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
				cellules[m_elementsARecevoir[voisin][i]]->recupereTamponPrimitives(m_tamponRec[0][voisin], count, eos, type);
      }
    } //Fin voisin
  }
}

//****************************************************************************
//******************** Methodes pour toutes les pentes ***********************
//****************************************************************************

void Parallel::initialiseCommunicationsPersistantesPentes()
{
	for (int voisin = 0; voisin < Ncpu; voisin++)	{
		if (m_estVoisin[voisin]) { //Si CPU voisin
			//Determination du nombre de variables a communiquer
			int nombreEnvoi = m_nombreVariablesPentes*m_nombreElementsAEnvoyerAVoisin[voisin];
			int nombreRecoi = m_nombreVariablesPentes*m_nombreElementsARecevoirDeVoisin[voisin];

			//Nouvelle requête d'envois et son buffer associe
			m_reqEnvoisPentes[0][voisin] = new MPI_Request;
			m_tamponEnvPentes[0][voisin] = new double[nombreEnvoi];
			MPI_Send_init(m_tamponEnvPentes[0][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisPentes[0][voisin]);

			//Nouvelle requête de receptions et son buffer associe
			m_reqReceptionsPentes[0][voisin] = new MPI_Request;
			m_tamponRecPentes[0][voisin] = new double[nombreRecoi];
			MPI_Recv_init(m_tamponRecPentes[0][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsPentes[0][voisin]);
		}
	}
}

//***********************************************************************

void Parallel::finaliseCommunicationsPersistantesPentes(const int &lvlMax)
{
	for (int lvl = 0; lvl <= lvlMax; lvl++) {
		for (int voisin = 0; voisin < Ncpu; voisin++)	{
			if (m_estVoisin[voisin]) { //Si CPU voisin
				MPI_Request_free(m_reqEnvoisPentes[lvl][voisin]);
				MPI_Request_free(m_reqReceptionsPentes[lvl][voisin]);
        delete m_reqEnvoisPentes[lvl][voisin];
        delete[] m_tamponEnvPentes[lvl][voisin];
        delete m_reqReceptionsPentes[lvl][voisin];
        delete[] m_tamponRecPentes[lvl][voisin];
			}
		}
		delete[] m_reqEnvoisPentes[lvl];
		delete[] m_tamponEnvPentes[lvl];
		delete[] m_reqReceptionsPentes[lvl];
		delete[] m_tamponRecPentes[lvl];
	}
  m_reqEnvoisPentes.clear();
  m_tamponEnvPentes.clear();
  m_reqReceptionsPentes.clear();
  m_tamponRecPentes.clear();
}

//***********************************************************************

void Parallel::communicationsPentes(Cellule **cellules)
{
	int count(0);
	MPI_Status status;

	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) {
			//Prepation des envois
			count = -1;
			for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
				cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponPentes(m_tamponEnvPentes[0][voisin], count);
			}
			
			//Requête d'envoi
			MPI_Start(m_reqEnvoisPentes[0][voisin]);
			//Requête de reception
			MPI_Start(m_reqReceptionsPentes[0][voisin]);
			//Attente
			MPI_Wait(m_reqEnvoisPentes[0][voisin], &status);
			MPI_Wait(m_reqReceptionsPentes[0][voisin], &status);

			//Receptions
			count = -1;
			for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
				cellules[m_elementsARecevoir[voisin][i]]->recupereTamponPentes(m_tamponRecPentes[0][voisin], count);
			}
		} //Fin voisin
	}
}

//****************************************************************************
//****************** Methodes pour une variable scalaire *********************
//****************************************************************************

void Parallel::initialiseCommunicationsPersistantesScalaire()
{
  int nombre;

  for (int voisin = 0; voisin < Ncpu; voisin++) {
    if (m_estVoisin[voisin]) { //Si CPU voisin
      //Determination du nombre de variables a communiquer
      nombre = 1; // 1 variable scalaire
      int nombreEnvoi = nombre*m_nombreElementsAEnvoyerAVoisin[voisin];
      int nombreRecoi = nombre*m_nombreElementsARecevoirDeVoisin[voisin];

      //Nouvelle requête d'envois et son buffer associe
      m_reqEnvoisScalaire[0][voisin] = new MPI_Request;
      m_tamponEnvScalaire[0][voisin] = new double[nombreEnvoi];
      MPI_Send_init(m_tamponEnvScalaire[0][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisScalaire[0][voisin]);

      //Nouvelle requête de receptions et son buffer associe
      m_reqReceptionsScalaire[0][voisin] = new MPI_Request;
      m_tamponRecScalaire[0][voisin] = new double[nombreRecoi];
      MPI_Recv_init(m_tamponRecScalaire[0][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsScalaire[0][voisin]);
    }
  }
}

//***********************************************************************

void Parallel::finaliseCommunicationsPersistantesScalaire(const int &lvlMax)
{
	for (int lvl = 0; lvl <= lvlMax; lvl++) {
		for (int voisin = 0; voisin < Ncpu; voisin++)	{
			if (m_estVoisin[voisin]) { //Si CPU voisin
				MPI_Request_free(m_reqEnvoisScalaire[lvl][voisin]);
				MPI_Request_free(m_reqReceptionsScalaire[lvl][voisin]);
        delete m_reqEnvoisScalaire[lvl][voisin];
        delete[] m_tamponEnvScalaire[lvl][voisin];
        delete m_reqReceptionsScalaire[lvl][voisin];
        delete[] m_tamponRecScalaire[lvl][voisin];
			}
		}
		delete[] m_reqEnvoisScalaire[lvl];
		delete[] m_tamponEnvScalaire[lvl];
		delete[] m_reqReceptionsScalaire[lvl];
		delete[] m_tamponRecScalaire[lvl];
	}
  m_reqEnvoisScalaire.clear();
  m_tamponEnvScalaire.clear();
  m_reqReceptionsScalaire.clear();
  m_tamponRecScalaire.clear();
}

//***********************************************************************

void Parallel::communicationsScalaire(Cellule **cellules, string nomScalaire)
{
  int count(0);
  MPI_Status status;

  for (int voisin = 0; voisin < Ncpu; voisin++) {
    if (m_estVoisin[voisin]) {
      //Prepation des envois
      count = -1;
      for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
        //Remplissage m_tamponEnvScalaire automatique selon nom de la variable
				cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponScalaire(m_tamponEnvScalaire[0][voisin], count, nomScalaire);
      }

      //Requête d'envoi
      MPI_Start(m_reqEnvoisScalaire[0][voisin]);
      //Requête de reception
      MPI_Start(m_reqReceptionsScalaire[0][voisin]);
      //Attente
      MPI_Wait(m_reqEnvoisScalaire[0][voisin], &status);
      MPI_Wait(m_reqReceptionsScalaire[0][voisin], &status);

      //Receptions
      count = -1;
      for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
        //Remplissage m_tamponRecGrad automatique selon coordonnees du gradient
				cellules[m_elementsARecevoir[voisin][i]]->recupereTamponScalaire(m_tamponRecScalaire[0][voisin], count, nomScalaire);
      }
    } //Fin voisin
  }
}

//****************************************************************************
//********************* Methodes pour les vecteurs ***************************
//****************************************************************************

void Parallel::initialiseCommunicationsPersistantesVecteur(const int &dim)
{
	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) { //Si CPU voisin
			//Determination du nombre de variables a communiquer, autant de variables que la dimension (1,2 ou 3)
			int nombreEnvoi = dim*m_nombreElementsAEnvoyerAVoisin[voisin];
			int nombreRecoi = dim*m_nombreElementsARecevoirDeVoisin[voisin];

			//Nouvelle requête d'envois et son buffer associe
			m_reqEnvoisVecteur[0][voisin] = new MPI_Request;
			m_tamponEnvVecteur[0][voisin] = new double[nombreEnvoi];
			MPI_Send_init(m_tamponEnvVecteur[0][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisVecteur[0][voisin]);

			//Nouvelle requête de receptions et son buffer associe
			m_reqReceptionsVecteur[0][voisin] = new MPI_Request;
			m_tamponRecVecteur[0][voisin] = new double[nombreRecoi];
			MPI_Recv_init(m_tamponRecVecteur[0][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsVecteur[0][voisin]);
		}
	}
}

//***********************************************************************

void Parallel::finaliseCommunicationsPersistantesVecteur(const int &lvlMax)
{
	for (int lvl = 0; lvl <= lvlMax; lvl++) {
		for (int voisin = 0; voisin < Ncpu; voisin++)	{
			if (m_estVoisin[voisin]) { //Si CPU voisin
				MPI_Request_free(m_reqEnvoisVecteur[lvl][voisin]);
				MPI_Request_free(m_reqReceptionsVecteur[lvl][voisin]);
        delete m_reqEnvoisVecteur[lvl][voisin];
        delete[] m_tamponEnvVecteur[lvl][voisin];
        delete m_reqReceptionsVecteur[lvl][voisin];
        delete[] m_tamponRecVecteur[lvl][voisin];
			}
		}
		delete[] m_reqEnvoisVecteur[lvl];
		delete[] m_tamponEnvVecteur[lvl];
		delete[] m_reqReceptionsVecteur[lvl];
		delete[] m_tamponRecVecteur[lvl];
	}
  m_reqEnvoisVecteur.clear();
  m_tamponEnvVecteur.clear();
  m_reqReceptionsVecteur.clear();
  m_tamponRecVecteur.clear();
}

//***********************************************************************

void Parallel::communicationsVecteur(Cellule **cellules, string nomVecteur, const int &dim, int num, int indice)
{
	int count(0);
	MPI_Status status;

	for (int voisin = 0; voisin < Ncpu; voisin++)	{
		if (m_estVoisin[voisin]) {
			//Prepation des envois
			count = -1;
			for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++)	{
			  //Remplissage m_tamponEnvVecteur automatique selon coordonnees du gradient
				cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponVecteur(m_tamponEnvVecteur[0][voisin], count, dim, nomVecteur, num, indice);
			}

			//Requête d'envoi
			MPI_Start(m_reqEnvoisVecteur[0][voisin]);
			//Requête de reception
			MPI_Start(m_reqReceptionsVecteur[0][voisin]);
			//Attente
			MPI_Wait(m_reqEnvoisVecteur[0][voisin], &status);
			MPI_Wait(m_reqReceptionsVecteur[0][voisin], &status);

			//Receptions
			count = -1;
			for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++)	{
			  //Remplissage m_tamponRecVecteur automatique selon coordonnees du gradient
				cellules[m_elementsARecevoir[voisin][i]]->recupereTamponVecteur(m_tamponRecVecteur[0][voisin], count, dim, nomVecteur, num, indice);
			}
		} //Fin voisin
	}
}

//****************************************************************************
//************ Methodes pour toutes les variables transportees ***************
//****************************************************************************

void Parallel::initialiseCommunicationsPersistantesTransports()
{
  for (int voisin = 0; voisin < Ncpu; voisin++) {
    if (m_estVoisin[voisin]) { //Si CPU voisin
      //Determination du nombre de variables a communiquer
      int nombreEnvoi = m_nombreVariablesTransports*m_nombreElementsAEnvoyerAVoisin[voisin];
      int nombreRecoi = m_nombreVariablesTransports*m_nombreElementsARecevoirDeVoisin[voisin];

      //Nouvelle requête d'envois et son buffer associe
      m_reqEnvoisTransports[0][voisin] = new MPI_Request;
      m_tamponEnvTransports[0][voisin] = new double[nombreEnvoi];
      MPI_Send_init(m_tamponEnvTransports[0][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisTransports[0][voisin]);

      //Nouvelle requête de receptions et son buffer associe
      m_reqReceptionsTransports[0][voisin] = new MPI_Request;
      m_tamponRecTransports[0][voisin] = new double[nombreRecoi];
      MPI_Recv_init(m_tamponRecTransports[0][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsTransports[0][voisin]);
    }
  }
}

//***********************************************************************

void Parallel::finaliseCommunicationsPersistantesTransports(const int &lvlMax)
{
  for (int lvl = 0; lvl <= lvlMax; lvl++) {
    for (int voisin = 0; voisin < Ncpu; voisin++) {
      if (m_estVoisin[voisin]) { //Si CPU voisin
        MPI_Request_free(m_reqEnvoisTransports[lvl][voisin]);
        MPI_Request_free(m_reqReceptionsTransports[lvl][voisin]);
        delete m_reqEnvoisTransports[lvl][voisin];
        delete[] m_tamponEnvTransports[lvl][voisin];
        delete m_reqReceptionsTransports[lvl][voisin];
        delete[] m_tamponRecTransports[lvl][voisin];
      }
    }
    delete[] m_reqEnvoisTransports[lvl];
    delete[] m_tamponEnvTransports[lvl];
    delete[] m_reqReceptionsTransports[lvl];
    delete[] m_tamponRecTransports[lvl];
  }
  m_reqEnvoisTransports.clear();
  m_tamponEnvTransports.clear();
  m_reqReceptionsTransports.clear();
  m_tamponRecTransports.clear();

}

//***********************************************************************

void Parallel::communicationsTransports(Cellule **cellules)
{
  int count(0);
  MPI_Status status;

  for (int voisin = 0; voisin < Ncpu; voisin++) {
    if (m_estVoisin[voisin]) {
      //Prepation des envois
      count = -1;
      for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
        cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponTransports(m_tamponEnvTransports[0][voisin], count);
      }

      //Requête d'envoi
      MPI_Start(m_reqEnvoisTransports[0][voisin]);
      //Requête de reception
      MPI_Start(m_reqReceptionsTransports[0][voisin]);
      //Attente
      MPI_Wait(m_reqEnvoisTransports[0][voisin], &status);
      MPI_Wait(m_reqReceptionsTransports[0][voisin], &status);

      //Receptions
      count = -1;
      for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
        cellules[m_elementsARecevoir[voisin][i]]->recupereTamponTransports(m_tamponRecTransports[0][voisin], count);
      }
    } //Fin voisin
  }
}

//****************************************************************************
//******************** Methodes pour les variables AMR ***********************
//****************************************************************************

void Parallel::initialiseCommunicationsPersistantesAMR(const int &nombreVariablesPrimitives, const int &nombreVariablesPentes, const int &nombreVariablesTransports, const int &dim, const int &lvlMax)
{
	if (Ncpu > 1) {
		m_nombreVariablesPrimitives = nombreVariablesPrimitives;
		m_nombreVariablesPentes = nombreVariablesPentes;
    m_nombreVariablesTransports = nombreVariablesTransports;
		//Initialisation des communications des variables primitives du modele resolu
		Calcul_Parallele.initialiseCommunicationsPersistantesPrimitives();
		//Initialisation des communications des pentes pour l'ordre 2
		Calcul_Parallele.initialiseCommunicationsPersistantesPentes();
		//Initialisation des communications necessaires aux physiques additionelles (vecteurs de dim=3)
		Calcul_Parallele.initialiseCommunicationsPersistantesVecteur(dim);
    //Initialisation des communications des variables transportees
    Calcul_Parallele.initialiseCommunicationsPersistantesTransports();
		//Initialisation des communications pour les variables AMR
		Calcul_Parallele.initialiseCommunicationsPersistantesXi();
		Calcul_Parallele.initialiseCommunicationsPersistantesSplit();
		Calcul_Parallele.initialiseCommunicationsPersistantesNombreCellulesGhost();
		//Initialisation des communications pour les niveaux superieurs a 0
		Calcul_Parallele.initialiseCommunicationsPersistantesNiveauxAMR(lvlMax);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::initialiseCommunicationsPersistantesNiveauxAMR(const int &lvlMax)
{
	//Extension des variables paralleles au niveau AMR maximum. On commence a 1, le niveau 0 etant deja initialise
	for (int lvl = 1; lvl <= lvlMax; lvl++) {
		m_tamponEnv.push_back(new double*[Ncpu]);
		m_tamponRec.push_back(new double*[Ncpu]);
		m_tamponEnvPentes.push_back(new double*[Ncpu]);
		m_tamponRecPentes.push_back(new double*[Ncpu]);
		m_tamponEnvVecteur.push_back(new double*[Ncpu]);
		m_tamponRecVecteur.push_back(new double*[Ncpu]);
    m_tamponEnvTransports.push_back(new double*[Ncpu]);
    m_tamponRecTransports.push_back(new double*[Ncpu]);
		m_tamponEnvXi.push_back(new double*[Ncpu]);
		m_tamponRecXi.push_back(new double*[Ncpu]);
		m_tamponEnvSplit.push_back(new bool*[Ncpu]);
		m_tamponRecSplit.push_back(new bool*[Ncpu]);

		m_reqEnvois.push_back(new MPI_Request*[Ncpu]);
		m_reqReceptions.push_back(new MPI_Request*[Ncpu]);
		m_reqEnvoisPentes.push_back(new MPI_Request*[Ncpu]);
		m_reqReceptionsPentes.push_back(new MPI_Request*[Ncpu]);
		m_reqEnvoisVecteur.push_back(new MPI_Request*[Ncpu]);
		m_reqReceptionsVecteur.push_back(new MPI_Request*[Ncpu]);
    m_reqEnvoisTransports.push_back(new MPI_Request*[Ncpu]);
    m_reqReceptionsTransports.push_back(new MPI_Request*[Ncpu]);
		m_reqEnvoisXi.push_back(new MPI_Request*[Ncpu]);
		m_reqReceptionsXi.push_back(new MPI_Request*[Ncpu]);
		m_reqEnvoisSplit.push_back(new MPI_Request*[Ncpu]);
		m_reqReceptionsSplit.push_back(new MPI_Request*[Ncpu]);

		for (int i = 0; i < Ncpu; i++) {
			m_tamponEnv[lvl][i] = NULL;
			m_tamponRec[lvl][i] = NULL;
			m_reqEnvois[lvl][i] = NULL;
			m_reqReceptions[lvl][i] = NULL;
			m_tamponEnvPentes[lvl][i] = NULL;
			m_tamponRecPentes[lvl][i] = NULL;
			m_reqEnvoisPentes[lvl][i] = NULL;
			m_reqReceptionsPentes[lvl][i] = NULL;
			m_tamponEnvVecteur[lvl][i] = NULL;
			m_tamponRecVecteur[lvl][i] = NULL;
			m_reqEnvoisVecteur[lvl][i] = NULL;
			m_reqReceptionsVecteur[lvl][i] = NULL;
      m_tamponEnvTransports[lvl][i] = NULL;
      m_tamponRecTransports[lvl][i] = NULL;
      m_reqEnvoisTransports[lvl][i] = NULL;
      m_reqReceptionsTransports[lvl][i] = NULL;
			m_tamponEnvXi[lvl][i] = NULL;
			m_tamponRecXi[lvl][i] = NULL;
			m_reqEnvoisXi[lvl][i] = NULL;
			m_reqReceptionsXi[lvl][i] = NULL;
			m_tamponEnvSplit[lvl][i] = NULL;
			m_tamponRecSplit[lvl][i] = NULL;
			m_reqEnvoisSplit[lvl][i] = NULL;
			m_reqReceptionsSplit[lvl][i] = NULL;
		}
	}

	//Initialisation des envois et receptions pour les couples de CPU voisins et pour chaque niveau AMR
	int nombreEnvoi(0);
	int nombreRecoi(0);

	for (int lvl = 1; lvl <= lvlMax; lvl++) {
		for (int voisin = 0; voisin < Ncpu; voisin++) {
			if (m_estVoisin[voisin]) { //Si CPU voisin
				//Variable Primitives
				//-------------------
				//Nouvelle requête d'envois et son buffer associe
				m_reqEnvois[lvl][voisin] = new MPI_Request;
				m_tamponEnv[lvl][voisin] = new double[nombreEnvoi];
				MPI_Send_init(m_tamponEnv[lvl][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvois[lvl][voisin]);

				//Nouvelle requête de receptions et son buffer associe
				m_reqReceptions[lvl][voisin] = new MPI_Request;
				m_tamponRec[lvl][voisin] = new double[nombreRecoi];
				MPI_Recv_init(m_tamponRec[lvl][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptions[lvl][voisin]);

				//Variable Pentes
				//---------------
				//Nouvelle requête d'envois et son buffer associe
				m_reqEnvoisPentes[lvl][voisin] = new MPI_Request;
				m_tamponEnvPentes[lvl][voisin] = new double[nombreEnvoi];
				MPI_Send_init(m_tamponEnvPentes[lvl][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisPentes[lvl][voisin]);

				//Nouvelle requête de receptions et son buffer associe
				m_reqReceptionsPentes[lvl][voisin] = new MPI_Request;
				m_tamponRecPentes[lvl][voisin] = new double[nombreRecoi];
				MPI_Recv_init(m_tamponRecPentes[lvl][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsPentes[lvl][voisin]);

				//Variable Vecteur
				//----------------
				//Nouvelle requête d'envois et son buffer associe
				m_reqEnvoisVecteur[lvl][voisin] = new MPI_Request;
				m_tamponEnvVecteur[lvl][voisin] = new double[nombreEnvoi];
				MPI_Send_init(m_tamponEnvVecteur[lvl][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisVecteur[lvl][voisin]);

				//Nouvelle requête de receptions et son buffer associe
				m_reqReceptionsVecteur[lvl][voisin] = new MPI_Request;
				m_tamponRecVecteur[lvl][voisin] = new double[nombreRecoi];
				MPI_Recv_init(m_tamponRecVecteur[lvl][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsVecteur[lvl][voisin]);

        //Variable Transportees
        //---------------------
        //Nouvelle requête d'envois et son buffer associe
        m_reqEnvoisTransports[lvl][voisin] = new MPI_Request;
        m_tamponEnvTransports[lvl][voisin] = new double[nombreEnvoi];
        MPI_Send_init(m_tamponEnvTransports[lvl][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisTransports[lvl][voisin]);

        //Nouvelle requête de receptions et son buffer associe
        m_reqReceptionsTransports[lvl][voisin] = new MPI_Request;
        m_tamponRecTransports[lvl][voisin] = new double[nombreRecoi];
        MPI_Recv_init(m_tamponRecTransports[lvl][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsTransports[lvl][voisin]);

				//Variable Xi
				//-----------
				//Nouvelle requête d'envois et son buffer associe
				m_reqEnvoisXi[lvl][voisin] = new MPI_Request;
				m_tamponEnvXi[lvl][voisin] = new double[nombreEnvoi];
				MPI_Send_init(m_tamponEnvXi[lvl][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisXi[lvl][voisin]);

				//Nouvelle requête de receptions et son buffer associe
				m_reqReceptionsXi[lvl][voisin] = new MPI_Request;
				m_tamponRecXi[lvl][voisin] = new double[nombreRecoi];
				MPI_Recv_init(m_tamponRecXi[lvl][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsXi[lvl][voisin]);

				//Variable Split
				//--------------
				//Nouvelle requête d'envois et son buffer associe
				m_reqEnvoisSplit[lvl][voisin] = new MPI_Request;
				m_tamponEnvSplit[lvl][voisin] = new bool[nombreEnvoi];
				MPI_Send_init(m_tamponEnvSplit[lvl][voisin], nombreEnvoi, MPI_C_BOOL, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisSplit[lvl][voisin]);

				//Nouvelle requête de receptions et son buffer associe
				m_reqReceptionsSplit[lvl][voisin] = new MPI_Request;
				m_tamponRecSplit[lvl][voisin] = new bool[nombreRecoi];
				MPI_Recv_init(m_tamponRecSplit[lvl][voisin], nombreRecoi, MPI_C_BOOL, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsSplit[lvl][voisin]);
			}
		}
	}
}

//***********************************************************************

void Parallel::miseAJourCommunicationsPersistantesLvl(int lvl, const int &dim)
{
	int nombreEnvoi(0), nombreRecoi(0);
	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) { //Si CPU voisin
			//---------------------------------------------------------------------
			//On vide dans un premier temps les variables d'envois et de receptions
			//---------------------------------------------------------------------
			MPI_Request_free(m_reqEnvois[lvl][voisin]);
			MPI_Request_free(m_reqReceptions[lvl][voisin]);
			MPI_Request_free(m_reqEnvoisPentes[lvl][voisin]);
			MPI_Request_free(m_reqReceptionsPentes[lvl][voisin]);
			MPI_Request_free(m_reqEnvoisVecteur[lvl][voisin]);
			MPI_Request_free(m_reqReceptionsVecteur[lvl][voisin]);
      MPI_Request_free(m_reqEnvoisTransports[lvl][voisin]);
      MPI_Request_free(m_reqReceptionsTransports[lvl][voisin]);
			MPI_Request_free(m_reqEnvoisXi[lvl][voisin]);
			MPI_Request_free(m_reqReceptionsXi[lvl][voisin]);
			MPI_Request_free(m_reqEnvoisSplit[lvl][voisin]);
			MPI_Request_free(m_reqReceptionsSplit[lvl][voisin]);

			delete m_reqEnvois[lvl][voisin];
			delete m_reqReceptions[lvl][voisin];
			delete m_reqEnvoisPentes[lvl][voisin];
			delete m_reqReceptionsPentes[lvl][voisin];
			delete m_reqEnvoisVecteur[lvl][voisin];
			delete m_reqReceptionsVecteur[lvl][voisin];
      delete m_reqEnvoisTransports[lvl][voisin];
      delete m_reqReceptionsTransports[lvl][voisin];
			delete m_reqEnvoisXi[lvl][voisin];
			delete m_reqReceptionsXi[lvl][voisin];
			delete m_reqEnvoisSplit[lvl][voisin];
			delete m_reqReceptionsSplit[lvl][voisin];

			delete[] m_tamponEnv[lvl][voisin];
			delete[] m_tamponRec[lvl][voisin];
			delete[] m_tamponEnvPentes[lvl][voisin];
			delete[] m_tamponRecPentes[lvl][voisin];
			delete[] m_tamponEnvVecteur[lvl][voisin];
			delete[] m_tamponRecVecteur[lvl][voisin];
      delete[] m_tamponEnvTransports[lvl][voisin];
      delete[] m_tamponRecTransports[lvl][voisin];
			delete[] m_tamponEnvXi[lvl][voisin];
			delete[] m_tamponRecXi[lvl][voisin];
			delete[] m_tamponEnvSplit[lvl][voisin];
			delete[] m_tamponRecSplit[lvl][voisin];

			//----------------------------------------------------------
			//On ecrit les nouvelles variables d'envois et de receptions
			//----------------------------------------------------------

			//Variable Primitives
			//-------------------
			nombreEnvoi = m_nombreVariablesPrimitives*m_tamponNombreElementsAEnvoyerAVoisin[voisin];
			nombreRecoi = m_nombreVariablesPrimitives*m_tamponNombreElementsARecevoirDeVoisin[voisin];
			//Nouvelle requête d'envois et son buffer associe
			m_reqEnvois[lvl][voisin] = new MPI_Request;
			m_tamponEnv[lvl][voisin] = new double[nombreEnvoi];
			MPI_Send_init(m_tamponEnv[lvl][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvois[lvl][voisin]);

			//Nouvelle requête de receptions et son buffer associe
			m_reqReceptions[lvl][voisin] = new MPI_Request;
			m_tamponRec[lvl][voisin] = new double[nombreRecoi];
			MPI_Recv_init(m_tamponRec[lvl][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptions[lvl][voisin]);

			//Variable Pentes
			//---------------
			nombreEnvoi = m_nombreVariablesPentes*m_tamponNombreElementsAEnvoyerAVoisin[voisin];
			nombreRecoi = m_nombreVariablesPentes*m_tamponNombreElementsARecevoirDeVoisin[voisin];
			//Nouvelle requête d'envois et son buffer associe
			m_reqEnvoisPentes[lvl][voisin] = new MPI_Request;
			m_tamponEnvPentes[lvl][voisin] = new double[nombreEnvoi];
			MPI_Send_init(m_tamponEnvPentes[lvl][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisPentes[lvl][voisin]);

			//Nouvelle requête de receptions et son buffer associe
			m_reqReceptionsPentes[lvl][voisin] = new MPI_Request;
			m_tamponRecPentes[lvl][voisin] = new double[nombreRecoi];
			MPI_Recv_init(m_tamponRecPentes[lvl][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsPentes[lvl][voisin]);

			//Variable Vecteur
			//----------------
			nombreEnvoi = dim*m_tamponNombreElementsAEnvoyerAVoisin[voisin];
			nombreRecoi = dim*m_tamponNombreElementsARecevoirDeVoisin[voisin];
			//Nouvelle requête d'envois et son buffer associe
			m_reqEnvoisVecteur[lvl][voisin] = new MPI_Request;
			m_tamponEnvVecteur[lvl][voisin] = new double[nombreEnvoi];
			MPI_Send_init(m_tamponEnvVecteur[lvl][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisVecteur[lvl][voisin]);

			//Nouvelle requête de receptions et son buffer associe
			m_reqReceptionsVecteur[lvl][voisin] = new MPI_Request;
			m_tamponRecVecteur[lvl][voisin] = new double[nombreRecoi];
			MPI_Recv_init(m_tamponRecVecteur[lvl][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsVecteur[lvl][voisin]);

      //Variable Transportees
      //---------------------
      nombreEnvoi = m_nombreVariablesTransports*m_tamponNombreElementsAEnvoyerAVoisin[voisin];
      nombreRecoi = m_nombreVariablesTransports*m_tamponNombreElementsARecevoirDeVoisin[voisin];
      //Nouvelle requête d'envois et son buffer associe
      m_reqEnvoisTransports[lvl][voisin] = new MPI_Request;
      m_tamponEnvTransports[lvl][voisin] = new double[nombreEnvoi];
      MPI_Send_init(m_tamponEnvTransports[lvl][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisTransports[lvl][voisin]);

      //Nouvelle requête de receptions et son buffer associe
      m_reqReceptionsTransports[lvl][voisin] = new MPI_Request;
      m_tamponRecTransports[lvl][voisin] = new double[nombreRecoi];
      MPI_Recv_init(m_tamponRecTransports[lvl][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsTransports[lvl][voisin]);

			//Variable Xi
			//-----------
			nombreEnvoi = m_tamponNombreElementsAEnvoyerAVoisin[voisin];
			nombreRecoi = m_tamponNombreElementsARecevoirDeVoisin[voisin];
			//Nouvelle requête d'envois et son buffer associe
			m_reqEnvoisXi[lvl][voisin] = new MPI_Request;
			m_tamponEnvXi[lvl][voisin] = new double[nombreEnvoi];
			MPI_Send_init(m_tamponEnvXi[lvl][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisXi[lvl][voisin]);

			//Nouvelle requête de receptions et son buffer associe
			m_reqReceptionsXi[lvl][voisin] = new MPI_Request;
			m_tamponRecXi[lvl][voisin] = new double[nombreRecoi];
			MPI_Recv_init(m_tamponRecXi[lvl][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsXi[lvl][voisin]);

			//Variable Split
			//--------------
			//Nouvelle requête d'envois et son buffer associe
			m_reqEnvoisSplit[lvl][voisin] = new MPI_Request;
			m_tamponEnvSplit[lvl][voisin] = new bool[nombreEnvoi];
			MPI_Send_init(m_tamponEnvSplit[lvl][voisin], nombreEnvoi, MPI_C_BOOL, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisSplit[lvl][voisin]);

			//Nouvelle requête de receptions et son buffer associe
			m_reqReceptionsSplit[lvl][voisin] = new MPI_Request;
			m_tamponRecSplit[lvl][voisin] = new bool[nombreRecoi];
			MPI_Recv_init(m_tamponRecSplit[lvl][voisin], nombreRecoi, MPI_C_BOOL, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsSplit[lvl][voisin]);
		}
	}
}

//***********************************************************************

void Parallel::finaliseAMR(const int &lvlMax)
{
	if (Ncpu > 1) {
		this->finaliseCommunicationsPersistantesPrimitives(lvlMax);
		this->finaliseCommunicationsPersistantesPentes(lvlMax);
		this->finaliseCommunicationsPersistantesVecteur(lvlMax);
    this->finaliseCommunicationsPersistantesTransports(lvlMax);
		this->finaliseCommunicationsPersistantesXi(lvlMax);
		this->finaliseCommunicationsPersistantesSplit(lvlMax);
		this->finaliseCommunicationsPersistantesNombreCellulesGhost();
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::initialiseCommunicationsPersistantesXi()
{
	int nombre(1);

	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) { //Si CPU voisin
			//Determination du nombre de variables a communiquer
			int nombreEnvoi = nombre*m_nombreElementsAEnvoyerAVoisin[voisin];
			int nombreRecoi = nombre*m_nombreElementsARecevoirDeVoisin[voisin];

			//Nouvelle requête d'envois et son buffer associe
			m_reqEnvoisXi[0][voisin] = new MPI_Request;
			m_tamponEnvXi[0][voisin] = new double[nombreEnvoi];
			MPI_Send_init(m_tamponEnvXi[0][voisin], nombreEnvoi, MPI_DOUBLE_PRECISION, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisXi[0][voisin]);

			//Nouvelle requête de receptions et son buffer associe
			m_reqReceptionsXi[0][voisin] = new MPI_Request;
			m_tamponRecXi[0][voisin] = new double[nombreRecoi];
			MPI_Recv_init(m_tamponRecXi[0][voisin], nombreRecoi, MPI_DOUBLE_PRECISION, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsXi[0][voisin]);
		}
	}
}

//***********************************************************************

void Parallel::finaliseCommunicationsPersistantesXi(const int &lvlMax)
{
	for (int lvl = 0; lvl <= lvlMax; lvl++) {
		for (int voisin = 0; voisin < Ncpu; voisin++) {
			if (m_estVoisin[voisin]) { //Si CPU voisin
				MPI_Request_free(m_reqEnvoisXi[lvl][voisin]);
				MPI_Request_free(m_reqReceptionsXi[lvl][voisin]);
        delete m_reqEnvoisXi[lvl][voisin];
        delete[] m_tamponEnvXi[lvl][voisin];
        delete m_reqReceptionsXi[lvl][voisin];
        delete[] m_tamponRecXi[lvl][voisin];
			}
		}
		delete[] m_reqEnvoisXi[lvl];
		delete[] m_tamponEnvXi[lvl];
		delete[] m_reqReceptionsXi[lvl];
		delete[] m_tamponRecXi[lvl];
	}
  m_reqEnvoisXi.clear();
  m_tamponEnvXi.clear();
  m_reqReceptionsXi.clear();
  m_tamponRecXi.clear();
}

//***********************************************************************

void Parallel::communicationsXi(Cellule **cellules, const int &lvl)
{
	int count(0);
	MPI_Status status;

	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) {
			//Prepation des envois
			count = -1;
			if (rang < voisin) { //Je suis CPU gauche
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					//Remplissage m_tamponEnvXi automatique
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponXiJeSuisCpuGauche(m_tamponEnvXi[lvl][voisin], count, lvl);
				}
			}
			else { //Je suis CPU droite
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					//Remplissage m_tamponEnvXi automatique
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponXiJeSuisCpuDroite(m_tamponEnvXi[lvl][voisin], count, lvl);
				}
			}

			//Requête d'envoi
			MPI_Start(m_reqEnvoisXi[lvl][voisin]);
			//Requête de reception
			MPI_Start(m_reqReceptionsXi[lvl][voisin]);
			//Attente
			MPI_Wait(m_reqEnvoisXi[lvl][voisin], &status);
			MPI_Wait(m_reqReceptionsXi[lvl][voisin], &status);

			//Receptions
			count = -1;
			for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
				//Remplissage m_tamponRecXi automatique
				cellules[m_elementsARecevoir[voisin][i]]->recupereTamponXi(m_tamponRecXi[lvl][voisin], count, lvl);
			}
		} //Fin voisin
	}
}

//***********************************************************************

void Parallel::initialiseCommunicationsPersistantesSplit()
{
	int nombre(1);

	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) { //Si CPU voisin
			//Determination du nombre de variables a communiquer
			int nombreEnvoi = nombre*m_nombreElementsAEnvoyerAVoisin[voisin];
			int nombreRecoi = nombre*m_nombreElementsARecevoirDeVoisin[voisin];

			//Nouvelle requête d'envois et son buffer associe
			m_reqEnvoisSplit[0][voisin] = new MPI_Request;
			m_tamponEnvSplit[0][voisin] = new bool[nombreEnvoi];
			MPI_Send_init(m_tamponEnvSplit[0][voisin], nombreEnvoi, MPI_C_BOOL, voisin, voisin, MPI_COMM_WORLD, m_reqEnvoisSplit[0][voisin]);

			//Nouvelle requête de receptions et son buffer associe
			m_reqReceptionsSplit[0][voisin] = new MPI_Request;
			m_tamponRecSplit[0][voisin] = new bool[nombreRecoi];
			MPI_Recv_init(m_tamponRecSplit[0][voisin], nombreRecoi, MPI_C_BOOL, voisin, rang, MPI_COMM_WORLD, m_reqReceptionsSplit[0][voisin]);
		}
	}
}

//***********************************************************************

void Parallel::finaliseCommunicationsPersistantesSplit(const int &lvlMax)
{
	for (int lvl = 0; lvl <= lvlMax; lvl++) {
		for (int voisin = 0; voisin < Ncpu; voisin++) {
			if (m_estVoisin[voisin]) { //Si CPU voisin
				MPI_Request_free(m_reqEnvoisSplit[lvl][voisin]);
				MPI_Request_free(m_reqReceptionsSplit[lvl][voisin]);
        delete m_reqEnvoisSplit[lvl][voisin];
        delete[] m_tamponEnvSplit[lvl][voisin];
        delete m_reqReceptionsSplit[lvl][voisin];
        delete[] m_tamponRecSplit[lvl][voisin];
			}
		}
		delete[] m_reqEnvoisSplit[lvl];
		delete[] m_tamponEnvSplit[lvl];
		delete[] m_reqReceptionsSplit[lvl];
		delete[] m_tamponRecSplit[lvl];
	}
  m_reqEnvoisSplit.clear();
  m_tamponEnvSplit.clear();
  m_reqReceptionsSplit.clear();
  m_tamponRecSplit.clear();
}

//***********************************************************************

void Parallel::communicationsSplit(Cellule **cellules, const int &lvl)
{
	int count(0);
	MPI_Status status;

	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) {
			//Prepation des envois
			count = -1;
			if (rang < voisin) { //Je suis CPU gauche
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					//Remplissage m_tamponEnvSplit automatique
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponSplitJeSuisCpuGauche(m_tamponEnvSplit[lvl][voisin], count, lvl);
				}
			}
			else { //Je suis CPU droite
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					//Remplissage m_tamponEnvSplit automatique
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponSplitJeSuisCpuDroite(m_tamponEnvSplit[lvl][voisin], count, lvl);
				}
			}

			//Requête d'envoi
			MPI_Start(m_reqEnvoisSplit[lvl][voisin]);
			//Requête de reception
			MPI_Start(m_reqReceptionsSplit[lvl][voisin]);
			//Attente
			MPI_Wait(m_reqEnvoisSplit[lvl][voisin], &status);
			MPI_Wait(m_reqReceptionsSplit[lvl][voisin], &status);

			//Receptions
			count = -1;
			for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
				//Remplissage m_tamponRecSplit automatique
				cellules[m_elementsARecevoir[voisin][i]]->recupereTamponSplit(m_tamponRecSplit[lvl][voisin], count, lvl);
			}
		} //Fin voisin
	}
}

//***********************************************************************

void Parallel::initialiseCommunicationsPersistantesNombreCellulesGhost()
{
	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) { //Si CPU voisin
			//Determination du nombre de variables a communiquer
			int nombreEnvoi = 1;
			int nombreRecoi = 1;

			//Nouvelle requête d'envois et son buffer associe
			m_reqNombreElementsAEnvoyerAVoisin[voisin] = new MPI_Request;
			m_tamponNombreElementsAEnvoyerAVoisin[voisin] = 0;
			MPI_Send_init(&m_tamponNombreElementsAEnvoyerAVoisin[voisin], nombreEnvoi, MPI_INT, voisin, voisin, MPI_COMM_WORLD, m_reqNombreElementsAEnvoyerAVoisin[voisin]);

			//Nouvelle requête de receptions et son buffer associe
			m_reqNombreElementsARecevoirDeVoisin[voisin] = new MPI_Request;
			m_tamponNombreElementsARecevoirDeVoisin[voisin] = 0;
			MPI_Recv_init(&m_tamponNombreElementsARecevoirDeVoisin[voisin], nombreRecoi, MPI_INT, voisin, rang, MPI_COMM_WORLD, m_reqNombreElementsARecevoirDeVoisin[voisin]);
		}
	}
}

//***********************************************************************

void Parallel::finaliseCommunicationsPersistantesNombreCellulesGhost()
{
	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) { //Si CPU voisin
			MPI_Request_free(m_reqNombreElementsAEnvoyerAVoisin[voisin]);
			MPI_Request_free(m_reqNombreElementsARecevoirDeVoisin[voisin]);
			delete m_reqNombreElementsAEnvoyerAVoisin[voisin];
			delete m_reqNombreElementsARecevoirDeVoisin[voisin];
		}
	}
	delete[] m_reqNombreElementsAEnvoyerAVoisin;
	delete[] m_reqNombreElementsARecevoirDeVoisin;
	delete[] m_tamponNombreElementsAEnvoyerAVoisin;
	delete[] m_tamponNombreElementsARecevoirDeVoisin;
}

//***********************************************************************

void Parallel::communicationsNombreCellulesGhost(Cellule **cellules, const int &lvl)
{
	MPI_Status status;

	for (int voisin = 0; voisin < Ncpu; voisin++)
	{
		if (m_estVoisin[voisin]) {
			//Prepation de l'envoi
			m_tamponNombreElementsAEnvoyerAVoisin[voisin] = 0;
			if (rang < voisin) { //Je suis CPU gauche
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					//Remplissage m_tamponEnvSplit automatique
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliNombreElementsAEnvoyerAVoisinJeSuisCpuGauche(m_tamponNombreElementsAEnvoyerAVoisin[voisin], lvl);
				}
			}
			else { //Je suis CPU droite
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					//Remplissage m_tamponEnvSplit automatique
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliNombreElementsAEnvoyerAVoisinJeSuisCpuDroite(m_tamponNombreElementsAEnvoyerAVoisin[voisin], lvl);
				}
			}

			//Requête d'envoi
			MPI_Start(m_reqNombreElementsAEnvoyerAVoisin[voisin]);
			//Requête de reception
			MPI_Start(m_reqNombreElementsARecevoirDeVoisin[voisin]);
			//Attente
			MPI_Wait(m_reqNombreElementsAEnvoyerAVoisin[voisin], &status);
			MPI_Wait(m_reqNombreElementsARecevoirDeVoisin[voisin], &status);

			//Pas de reception supplementaire a gerer

		} //Fin voisin
	}
}

//***********************************************************************

void Parallel::communicationsPrimitivesAMR(Cellule **cellules, Eos **eos, const int &lvl, Prim type)
{
	int count(0);
	MPI_Status status;

	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) {
			//Prepation des envois
			count = -1;
			if (rang < voisin) { //Je suis CPU gauche
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++)	{
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponPrimitivesAMRjeSuisCpuGauche(m_tamponEnv[lvl][voisin], count, lvl, type);
				}
			}
			else { //Je suis CPU droite
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++)	{
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponPrimitivesAMRjeSuisCpuDroite(m_tamponEnv[lvl][voisin], count, lvl, type);
				}
			}

			//Requête d'envoi
			MPI_Start(m_reqEnvois[lvl][voisin]);
			//Requête de reception
			MPI_Start(m_reqReceptions[lvl][voisin]);
			//Attente
			MPI_Wait(m_reqEnvois[lvl][voisin], &status);
			MPI_Wait(m_reqReceptions[lvl][voisin], &status);

			//Receptions
			count = -1;
			for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
				cellules[m_elementsARecevoir[voisin][i]]->recupereTamponPrimitivesAMR(m_tamponRec[lvl][voisin], count, lvl, eos, type);
			}
		} //Fin voisin
	}
}

//***********************************************************************

void Parallel::communicationsPentesAMR(Cellule **cellules, const int &lvl)
{
	int count(0);
	MPI_Status status;

	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) {
			//Prepation des envois
			count = -1;
			if (rang < voisin) { //Je suis CPU gauche
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponPentesAMRjeSuisCpuGauche(m_tamponEnvPentes[lvl][voisin], count, lvl);
				}
			}
			else { //Je suis CPU droite
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponPentesAMRjeSuisCpuDroite(m_tamponEnvPentes[lvl][voisin], count, lvl);
				}
			}

			//Requête d'envoi
			MPI_Start(m_reqEnvoisPentes[lvl][voisin]);
			//Requête de reception
			MPI_Start(m_reqReceptionsPentes[lvl][voisin]);
			//Attente
			MPI_Wait(m_reqEnvoisPentes[lvl][voisin], &status);
			MPI_Wait(m_reqReceptionsPentes[lvl][voisin], &status);

			//Receptions
			count = -1;
			for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
				cellules[m_elementsARecevoir[voisin][i]]->recupereTamponPentesAMR(m_tamponRecPentes[lvl][voisin], count, lvl);
			}
		} //Fin voisin
	}
}

//***********************************************************************

void Parallel::communicationsScalaireAMR(Cellule **cellules, string nomScalaire, const int &lvl)
{
	int count(0);
	MPI_Status status;

	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) {
			//Prepation des envois
			count = -1;
			if (rang < voisin) { //Je suis CPU gauche
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					//Remplissage m_tamponEnvScalaire automatique selon nom de la variable
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponScalaireAMRjeSuisCpuGauche(m_tamponEnvScalaire[lvl][voisin], count, lvl, nomScalaire);
				}
			}
			else { //Je suis CPU droite
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					//Remplissage m_tamponEnvScalaire automatique selon nom de la variable
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponScalaireAMRjeSuisCpuDroite(m_tamponEnvScalaire[lvl][voisin], count, lvl, nomScalaire);
				}
			}

			//Requête d'envoi
			MPI_Start(m_reqEnvoisScalaire[lvl][voisin]);
			//Requête de reception
			MPI_Start(m_reqReceptionsScalaire[lvl][voisin]);
			//Attente
			MPI_Wait(m_reqEnvoisScalaire[lvl][voisin], &status);
			MPI_Wait(m_reqReceptionsScalaire[lvl][voisin], &status);

			//Receptions
			count = -1;
			for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
				//Remplissage m_tamponRecGrad automatique selon coordonnees du gradient
				cellules[m_elementsARecevoir[voisin][i]]->recupereTamponScalaireAMR(m_tamponRecScalaire[lvl][voisin], count, lvl, nomScalaire);
			}
		} //Fin voisin
	}
}

//***********************************************************************

void Parallel::communicationsVecteurAMR(Cellule **cellules, string nomVecteur, const int &dim, const int &lvl, int num, int indice)
{
	int count(0);
	MPI_Status status;

	for (int voisin = 0; voisin < Ncpu; voisin++) {
		if (m_estVoisin[voisin]) {
			//Prepation des envois
			count = -1;
			if (rang < voisin) { //Je suis CPU gauche
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					//Remplissage m_tamponEnvVecteur automatique selon coordonnees du gradient
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponVecteurAMRjeSuisCpuGauche(m_tamponEnvVecteur[lvl][voisin], count, lvl, dim, nomVecteur, num, indice);
				}
			}
			else { //Je suis CPU droite
				for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
					//Remplissage m_tamponEnvVecteur automatique selon coordonnees du gradient
					cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponVecteurAMRjeSuisCpuDroite(m_tamponEnvVecteur[lvl][voisin], count, lvl, dim, nomVecteur, num, indice);
				}
			}

			//Requête d'envoi
			MPI_Start(m_reqEnvoisVecteur[lvl][voisin]);
			//Requête de reception
			MPI_Start(m_reqReceptionsVecteur[lvl][voisin]);
			//Attente
			MPI_Wait(m_reqEnvoisVecteur[lvl][voisin], &status);
			MPI_Wait(m_reqReceptionsVecteur[lvl][voisin], &status);

			//Receptions
			count = -1;
			for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
				//Remplissage m_tamponRecVecteur automatique selon coordonnees du gradient
				cellules[m_elementsARecevoir[voisin][i]]->recupereTamponVecteurAMR(m_tamponRecVecteur[lvl][voisin], count, lvl, dim, nomVecteur, num, indice);
			}
		} //Fin voisin
	}
}

//***********************************************************************

void Parallel::communicationsTransportsAMR(Cellule **cellules, const int &lvl)
{
  int count(0);
  MPI_Status status;

  for (int voisin = 0; voisin < Ncpu; voisin++) {
    if (m_estVoisin[voisin]) {
      //Prepation des envois
      count = -1;
      if (rang < voisin) { //Je suis CPU gauche
        for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
          cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponTransportsAMRjeSuisCpuGauche(m_tamponEnvTransports[lvl][voisin], count, lvl);
        }
      }
      else { //Je suis CPU droite
        for (int i = 0; i < m_nombreElementsAEnvoyerAVoisin[voisin]; i++) {
          cellules[m_elementsAEnvoyer[voisin][i]]->rempliTamponTransportsAMRjeSuisCpuDroite(m_tamponEnvTransports[lvl][voisin], count, lvl);
        }
      }

      //Requête d'envoi
      MPI_Start(m_reqEnvoisTransports[lvl][voisin]);
      //Requête de reception
      MPI_Start(m_reqReceptionsTransports[lvl][voisin]);
      //Attente
      MPI_Wait(m_reqEnvoisTransports[lvl][voisin], &status);
      MPI_Wait(m_reqReceptionsTransports[lvl][voisin], &status);

      //Receptions
      count = -1;
      for (int i = 0; i < m_nombreElementsARecevoirDeVoisin[voisin]; i++) {
        cellules[m_elementsARecevoir[voisin][i]]->recupereTamponTransportsAMR(m_tamponRecTransports[lvl][voisin], count, lvl);
      }
    } //Fin voisin
  }
}

//***********************************************************************