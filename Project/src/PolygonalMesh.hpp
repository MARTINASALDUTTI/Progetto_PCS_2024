#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

namespace PolygonalMeshLibrary{
struct PolygonalMesh{

    unsigned int Num0DsCell;
    //std::unordered_map<unsigned int, Eigen::Vector3d> coord0DsCellMap;
    //std::vector<unsigned int> Id0DsCell; secondo me non ha senso creare questo vector ma chiedi
    std::vector<Eigen::Vector3d> coord0DsCell;
    //std::unordered_map<Eigen::Vector3d*, unsigned int> coord0DsCellMap; //ha senso??

    /*Eigen::Vector3d è un tipo non copiabile per impostazione predefinita, il che significa che non puoi copiarlo
     * direttamente usando il costruttore di copia o l'operatore di assegnazione.
     * Una soluzione a questo problema potrebbe essere utilizzare un puntatore a Eigen::Vector3d anziché il vettore
     *  stesso come chiave della mappa. In questo modo, eviterai di dover copiare i vettori quando li inserisci nella mappa.
     */

    /*
    L'uso di std::find ha una complessità lineare O(n), il che significa che per vettori di grandi dimensioni, questo metodo
    potrebbe non essere molto efficiente. Se la tua applicazione richiede molte ricerche e inserimenti unici, considerare
    l'uso di un std::unordered_set<unsigned int> o un std::set<unsigned int>, che offrono inserimenti e ricerche più
    efficienti (in media O(1) per std::unordered_set e O(log n) per std::set).
    dobbiamo vedere differenza set and unordered_set
    */

    unsigned int Num1DsCell;
    //std::unordered_map<unsigned int, std::array<unsigned int, 2>> coord1DsCellMap;
    //std::vector<unsigned int> Id1DsCell;
    std::vector<std::array<unsigned int, 2>> coord1DsCell;

    unsigned int Num2DsCell;
    std::vector<std::vector<unsigned int>> Id2DsCell;
    std::vector<std::array<unsigned int, 3>> Cell2DssVertices = {}; ///< Cell2D Vertices indices, size 1 x NumberCell2DVertices[NumberCell2D]
    std::vector<std::array<unsigned int, 3>> Cell2DsEdges = {}; ///< Cell2D Cell1D indices, size 1 x NumberCell2DEdges[NumberCell2D]
};
}

