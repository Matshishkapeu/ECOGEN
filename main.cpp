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

#include "src/Run.h"
#include "src/Erreurs.h"
#include "src/libTierces/tinyxml2.h"

using namespace std;
using namespace tinyxml2;

void afficheEntete();

//***********************************************************************

int main(int argc, char* argv[])
{
  Run* run(0);

  // set_terminate(Global::arretApresErreur);

  //Initialisation parallele
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rang);
  MPI_Comm_size(MPI_COMM_WORLD, &Ncpu);

  if(rang==0) afficheEntete();
  MPI_Barrier(MPI_COMM_WORLD);

  //Parsing du fichier XML par la bibliotheque tinyxml2
  //---------------------------------------------------
  stringstream nomFichier("ECOGEN.xml");
  XMLDocument xmlEcogen;
  XMLError erreur(xmlEcogen.LoadFile(nomFichier.str().c_str())); //Le fichier est parse ici
  //if (erreur != XML_SUCCESS) throw ErreurXML(nomFichier.str(), __FILE__, __LINE__);
  //Recuperation racine du document XML
  XMLNode *xmlNode = xmlEcogen.FirstChildElement("ecogen");
  //if (xmlNode == NULL) throw ErreurXMLRacine("ecogen", nomFichier.str(), __FILE__, __LINE__);

  //Boucle sur les cas tests a executer
  //-----------------------------------
  int numCasTest(0);
  XMLElement *elementCasTest = xmlNode->FirstChildElement("casTest");
  while (elementCasTest != NULL) {
    try {
      XMLNode* xmlNode2 = elementCasTest->FirstChild();
      if (xmlNode2 == NULL) throw ErreurXMLElement("casTest", nomFichier.str(), __FILE__, __LINE__);
      XMLText* xmlText = xmlNode2->ToText();
      if (xmlText == NULL) throw ErreurXMLElement("casTest", nomFichier.str(), __FILE__, __LINE__);

      //1) Creation du cas Test
      numCasTest++;
      run = new Run(xmlText->Value(), numCasTest);
      MPI_Barrier(MPI_COMM_WORLD);
      if (rang == 0) {
        cout << "************************************************************" << endl;
        cout << "          EXECUTION DU CAS TEST NUMERO : " << numCasTest << endl;
        cout << "************************************************************" << endl;
        cout << "T" << numCasTest << " | Cas test : " << xmlText->Value() << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      //2) Execution du cas test
      //if (rang == 0) { cout << "attente lancement test" << endl;  system("pause"); }
      //MPI_Barrier(MPI_COMM_WORLD);
      run->initialisation(argc, argv);
      run->resolution();
      //if (rang == 0) { cout << "attente fin test" << endl; system("pause"); }
      //MPI_Barrier(MPI_COMM_WORLD);
      run->finalise();
      //3) Suppression du cas test
      delete run;
      //if (rang == 0) { cout << "test fini" << endl; system("pause"); }
      //MPI_Barrier(MPI_COMM_WORLD);
    }
    //Gestion Exceptions
    //------------------
    catch (ErreurXML &) {
      if (!run) {
        run->finalise();
        delete run;
      }
    }
    catch (ErreurECOGEN &e) {
      cerr << e.infoErreur() << endl;
      //cerr << "T" << numCasTest << " | " << e.infoErreur() << endl;
      if (!run) {
        run->finalise();
        delete run;
      }
    }
    elementCasTest = elementCasTest->NextSiblingElement("casTest");
  }//Fin boucle de cas tests
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}

//***********************************************************************

void afficheEntete()
{
  cout << "************************************************************" << endl;
  cout << "    ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. " << endl;
  cout << "    | .-'   .' .')   / .-. )  .' .'     | .-'    |  \\| | " << endl;
  cout << "    | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | " << endl;
  cout << "    | .-'   \\  \\     | | | |  \\  \\ ( _) | .-'    | |\\  | " << endl;
  cout << "    |  `--.  \\  `-.  \\ `-' /   \\  `-) ) |  `--.  | | |)| " << endl;
  cout << "    /( __.'   \\____\\  )---'    )\\____/  /( __.'  /(  (_) " << endl;
  cout << "   (__)              (_)      (__)     (__)     (__)     " << endl;
  cout << "************************************************************" << endl;
}

//***********************************************************************