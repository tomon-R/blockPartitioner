#ifndef BLOCK_PARTITIONER_HPP
#define BLOCK_PARTITIONER_HPP

#pragma once
#include <cstddef>
#include <vector>
#include <array>
#include <string>
#include <functional>
#include <map>
#include <tuple>

namespace partitioner {

// ───── Values ─────
struct XYZ           { double x, y, z;   };                 
struct Index3        { size_t i, j, k;   };            // Grid index

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
    * Finds a nodeId by idx. 
    * @return nodeId
    */
    Node* find(Index3 idx) const;
    const std::vector<std::vector<std::vector<Node*>>>& getNodes() const;
    void setNodes(std::vector<std::vector<std::vector<Node*>>>&& n) {
        nodes = std::move(n);
    }
    void setBoundaries(double minX, double maxX, double minY, double maxY, double minZ, double maxZ) {
        this->minX = minX;
        this->maxX = maxX;
        this->minY = minY;
        this->maxY = maxY;
        this->minZ = minZ;
        this->maxZ = maxZ;
    }
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
     * It works because it assumes mesh to be structured. It can't apply to free mesh.
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

    static std::tuple<bool, std::string> ifBlockSizesEqual(const std::vector<Subdomain>& subdomainAssets, const double tolerance);
    static std::tuple<bool, std::string> ifNodeSizesEqual(const std::vector<Subdomain>& subdomainAssets);
    static std::tuple<bool, std::string> ifNodeOrdersEqual(const std::vector<Subdomain>& subdomainAssets);
    static std::tuple<bool, std::string> ifNodePositionsEqual(const std::vector<Subdomain>& subdomainAssets);
    static std::tuple<bool, std::string> ifAllNodesPossessed(const Domain& domain);


public:
    static std::tuple<bool, std::string> ifBlocksValid(const Domain& domain, const double tolerance);
    static std::tuple<bool, std::string> canDivideEvenly(const Mesh& mesh, int nx, int ny, int nz);

};

} // namespace partitioner

class Args {
    std::vector<std::string> positionals;
    std::map<std::string, std::string> options;

public:
    Args parse(int argc, char* argv[], size_t requiredPositionals) ;
    std::string getOpt(const std::string& key) const;
    const std::vector<std::string>& getPos() const { 
        return positionals; 
    }
};


#endif // BLOCK_PARTITIONER_HPP