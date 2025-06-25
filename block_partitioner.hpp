#ifndef BLOCK_PARTITIONER_HPP
#define BLOCK_PARTITIONER_HPP

#pragma once
#include <cstddef>
#include <vector>
#include <array>
#include <string>
#include <functional>

namespace partitioner {

// ───── Values ─────
struct XYZ           { double x, y, z;   };                 
struct Index3        { size_t i, j, k;   };            // Grid index
struct Tolerance     { double tolerance; };

// ───── Entities & Aggregates ─────
struct Node {                                        
    int   id;
    XYZ  pos;
    Index3 globalIdx;
    int subId;
    Index3 localIdx;
};

class Mesh {                                          
    std::vector<std::vector<std::vector<Node*>>> nodes;

    double minX, maxX, minY, maxY, minZ, maxZ;

public:
    /*
    * Finds a nodeId by idx. Costs N(O1).
    * @return nodeId
    */
    Node* find(Index3 idx) const;
    const std::vector<std::vector<std::vector<Node*>>>& getNodes() const;
};

struct Subdomain {
    int                   id;
    size_t                size;
    std::vector<Node*>    nodes;
    Mesh                  localMesh;
};

/*
* Takes account of memories management.
*/
struct Domain {             
    std::vector<Node>          nodes; 
    Mesh                       globalMesh;             
    std::vector<Subdomain>     blocks;
};

// ───── Services ─────

class FileIoService {
public:
    /*
    * Loads a node.dat-format file and returns a Domain.
    */
    static std::vector<Node> load(const std::string path);
    /*
    * Writes out the part.dat-format file
    */
    static void write(const std::string path, const Domain& domain);
};

class MeshingService {
public:
    /**
     * Builds a 3D array with structured positioned node assets
     */
    static Mesh structuredMesh(const std::vector<Node>& nodes);
};

class PartitionService {
public:
    /*
    * Block partition
    */
    static std::vector<Subdomain> block(const Mesh& mesh, int nx, int ny, int nz);
};

class ValidationService {

    static bool ifBlockSizesEqual(const std::vector<Subdomain>& subdomainAssets, const Tolerance tolerance);
    static bool ifNodeSizesEqual(const std::vector<Subdomain>& subdomainAssets);
    static bool ifNodeOrdersEqual(const std::vector<Subdomain>& subdomainAssets);
    static bool ifNodePositionsEqual(const std::vector<Subdomain>& subdomainAssets);

public:
    static bool ifBlocksValid(const Domain& domain);

};


class Util {
public:
    template <typename T>
    static void qSort(std::vector<T>& array, std::function<bool(const T, const T)> judge) {
        std::sort(array.begin(), array.end(), judge);
    }
};

} // namespace partitioner

#endif // BLOCK_PARTITIONER_HPP