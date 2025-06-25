#include "block_partitioner.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: block_partitioner NX NY NZ [-ig input] [-d output] [-e eps]" << std::endl;
        return 1;
    }

    int nx = std::stoi(argv[1]);
    int ny = std::stoi(argv[2]);
    int nz = std::stoi(argv[3]);

    std::string input_path = "./node.dat";
    std::string output_dir = "./";
    double epsilon = 0.0;

    for (int i = 4; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-ig" && i + 1 < argc) {
            input_path = argv[++i];
        } else if (arg == "-d" && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (arg == "-e" && i + 1 < argc) {
            epsilon = std::stod(argv[++i]);
        }
    }

    std::cout << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
    std::cout << "input = " << input_path << ", output = " << output_dir << ", eps = " << epsilon << std::endl;


    using namespace partitioner;

    // Load domain from file (dummy)
    Domain domain;
    
    domain.nodes = FileIoService::load("node.dat");

    // Build mesh
    domain.globalMesh = MeshingService::structuredMesh(domain.nodes);

    // Partition
    domain.blocks = PartitionService::block(domain.globalMesh, nx, ny, nz);

    // Validate
    if (ValidationService::ifBlocksValid(domain)) {
        std::cout << "Partitioning is valid.\n";
    } else {
        std::cerr << "Partitioning validation failed.\n";
    }

    // Write output
    FileIoService::write("part.dat", domain);

    return 0;
}



namespace partitioner {

Node* Mesh::find(Index3 idx) const {
    return nodes[idx.i][idx.j][idx.k];
}

const std::vector<std::vector<std::vector<Node*>>>& Mesh::getNodes() const {
    return nodes;
}

Domain FileIoService::load(const std::string path) {
    // Dummy loader - replace with actual parsing
    Domain domain;
    // example node
    domain.nodes.push_back({0, {0.0, 0.0, 0.0}, {0, 0, 0}, 0, {0, 0, 0}});
    return domain;
}

void FileIoService::write(const std::string path, const Domain& domain) {
    std::ofstream file(path);
    for (const auto& node : domain.nodes) {
        file << node.id << ", " << node.pos.x << ", " << node.pos.y << ", " << node.pos.z << "\n";
    }
}

Mesh MeshingService::structuredMesh(const std::vector<Node>& nodes) {
    Mesh mesh;
    // Dummy mesh assignment - fill mesh.nodes properly in practice
    return mesh;
}

std::vector<Subdomain> PartitionService::block(const Mesh& mesh) {
    std::vector<Subdomain> result;
    // Dummy partitioning - divide into 2 for example
    Subdomain s1{0, 0, {}, mesh};
    Subdomain s2{1, 0, {}, mesh};
    result.push_back(s1);
    result.push_back(s2);
    return result;
}


// ValidationService private static methods
bool ValidationService::ifBlockSizesEqual(const std::vector<Subdomain>& subdomainAssets, const Tolerance tolerance) {
    if (subdomainAssets.empty()) return true;
    const size_t baseSize = subdomainAssets[0].size;
    for (const auto& sub : subdomainAssets) {
        if (std::abs(static_cast<double>(sub.size - baseSize)) > tolerance.tolerance * baseSize) {
            return false;
        }
    }
    return true;
}

bool ValidationService::ifNodeSizesEqual(const std::vector<Subdomain>& subdomainAssets) {
    if (subdomainAssets.empty()) return true;
    const size_t baseSize = subdomainAssets[0].nodes.size();
    for (const auto& sub : subdomainAssets) {
        if (sub.nodes.size() != baseSize) return false;
    }
    return true;
}

bool ValidationService::ifNodeOrdersEqual(const std::vector<Subdomain>& subdomainAssets) {
    if (subdomainAssets.empty()) return true;
    const auto& ref = subdomainAssets[0].nodes;
    for (const auto& sub : subdomainAssets) {
        for (size_t i = 0; i < ref.size(); ++i) {
            if (ref[i]->id != sub.nodes[i]->id) return false;
        }
    }
    return true;
}

bool ValidationService::ifNodePositionsEqual(const std::vector<Subdomain>& subdomainAssets) {
    if (subdomainAssets.empty()) return true;
    const auto& ref = subdomainAssets[0].nodes;
    for (const auto& sub : subdomainAssets) {
        for (size_t i = 0; i < ref.size(); ++i) {
            const auto& a = ref[i]->pos;
            const auto& b = sub.nodes[i]->pos;
            if (a.x != b.x || a.y != b.y || a.z != b.z) return false;
        }
    }
    return true;
}

bool ValidationService::ifBlocksValid(const Domain& domain) {
    const auto& blocks = domain.blocks;
    return ifBlockSizesEqual(blocks, {0.01}) &&
           ifNodeSizesEqual(blocks) &&
           ifNodeOrdersEqual(blocks) &&
           ifNodePositionsEqual(blocks);
}

} // namespace partitioner
