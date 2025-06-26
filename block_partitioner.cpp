#include "block_partitioner.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>


using namespace partitioner;

int main(int argc, char* argv[]) {
    Args args = Args().parse(argc, argv, 3);

    std::vector<std::string>& positionals = args.getPos();
    int nx = std::stoi(positionals[0]);
    int ny = std::stoi(positionals[1]);
    int nz = std::stoi(positionals[2]);

    std::string input_path = args.getOpt("-ig").empty() ? args.getOpt("-ig") : "./node.dat";
    std::string output_dir = args.getOpt("-d").empty() ? args.getOpt("-d") : "./part.dat";
    double epsilon = args.getOpt("-e").empty() ? std::stod(args.getOpt("-e")) : 0.0;

    // Load domain from file (dummy)
    Domain domain;
    
    domain.nodes = FileIoService::load(input_path);

    // Build mesh
    domain.globalMesh = MeshingService::structuredMesh(domain.nodes);

    // Partition
    domain.blocks = PartitionService::block(domain.globalMesh, nx, ny, nz);

    bool isValid;
    std::string errorMessage;
    std::tie(isValid, errorMessage) = ValidationService::ifBlocksValid(domain);

    // Validate
    if (isValid) {
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

std::vector<Node> FileIoService::load(const std::string path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + path);
    }

    size_t n;
    file >> n;

    std::vector<Node> nodes;
    nodes.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        double x, y, z;
        file >> x >> y >> z;
        Node node;
        node.id = static_cast<int>(i);
        node.pos = {x, y, z};
        node.globalIdx = {0, 0, 0};  // structuredMeshで上書きされる
        node.subId = -1;
        node.localIdx = {0, 0, 0};
        nodes.push_back(node);
    }

    return nodes;
}

void FileIoService::write(const std::string path, const Domain& domain) {
    std::ofstream file(path);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + path);
    }

    const auto& blocks = domain.blocks;
    file << blocks.size() << "\n";
    for (const auto& block : blocks) {
        file << block.id << " " << block.nodes.size();
        for (const auto* node : block.nodes) {
            file << " " << node->id;
        }
        file << "\n";
    }
}

Mesh MeshingService::structuredMesh(const std::vector<Node>& nodes) {
    double TOLERANCE = 1e-8;
    Mesh mesh;
    if (nodes.empty()) return mesh;

    // 1. x, y, z 座標を収集
    std::vector<double> xs, ys, zs;
    for (const auto& node : nodes) {
        xs.push_back(node.pos.x);
        ys.push_back(node.pos.y);
        zs.push_back(node.pos.z);
    }

    // 2. 重複を除いてソート
    auto unique_sorted = [](std::vector<double>& v, double tolerance) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end(), [](double a, double b) {
            return std::abs(a - b) < tolerance;
        }), v.end());
    };

    unique_sorted(xs, TOLERANCE);
    unique_sorted(ys, TOLERANCE);
    unique_sorted(zs, TOLERANCE);

    size_t nx = xs.size();
    size_t ny = ys.size();
    size_t nz = zs.size();

    // 3. メッシュ用配列を初期化
    std::vector<std::vector<std::vector<Node*>>> grid(
        nx, std::vector<std::vector<Node*>>(
            ny, std::vector<Node*>(nz, nullptr)
        )
    );

    // 4. 各ノードにグリッドインデックスを割り当てて配置
    for (auto& node : const_cast<std::vector<Node>&>(nodes)) {
        auto find_idx = [](const std::vector<double>& v, double val, double tolerance) -> size_t {
            for (size_t i = 0; i < v.size(); ++i) {
                if (std::abs(v[i] - val) < tolerance) return i;
            }
            throw std::runtime_error("Coordinate not found in grid axis");
        };

        size_t i = find_idx(xs, node.pos.x, TOLERANCE);
        size_t j = find_idx(ys, node.pos.y, TOLERANCE);
        size_t k = find_idx(zs, node.pos.z, TOLERANCE);

        node.globalIdx = {i, j, k};
        grid[i][j][k] = &node;
    }

    // 5. メッシュを構築
    mesh = Mesh();
    mesh.setNodes(std::move(grid));

    return mesh;
}

std::vector<Subdomain> PartitionService::block(const Mesh& mesh, int nx, int ny, int nz) {
    const auto& grid = mesh.getNodes();
    size_t Nx = grid.size();
    size_t Ny = grid[0].size();
    size_t Nz = grid[0][0].size();

    size_t dx = Nx / nx;
    size_t dy = Ny / ny;
    size_t dz = Nz / nz;

    std::vector<Subdomain> blocks;
    int blockId = 0;

    for (int sx = 0; sx < nx; ++sx) {
        for (int sy = 0; sy < ny; ++sy) {
            for (int sz = 0; sz < nz; ++sz) {
                Subdomain sub;
                sub.id = blockId++;
                sub.nodes.clear();

                for (size_t i = sx * dx; i < (sx + 1) * dx; ++i) {
                    for (size_t j = sy * dy; j < (sy + 1) * dy; ++j) {
                        for (size_t k = sz * dz; k < (sz + 1) * dz; ++k) {
                            Node* node = grid[i][j][k];
                            node->subId = sub.id;
                            node->localIdx = {i - sx * dx, j - sy * dy, k - sz * dz};
                            sub.nodes.push_back(node);
                        }
                    }
                }

                sub.size = sub.nodes.size();
                sub.localMesh = {};
                blocks.push_back(std::move(sub));
            }
        }
    }

    return blocks;
}


// ValidationService private static methods
std::tuple<bool, std::string> ValidationService::ifBlockSizesEqual(const std::vector<Subdomain>& subdomainAssets, const Tolerance tolerance) {
    if (subdomainAssets.empty()) return std::forward_as_tuple(true, "");
    const size_t baseSize = subdomainAssets[0].size;
    for (const auto& sub : subdomainAssets) {
        if (std::abs(static_cast<double>(sub.size - baseSize)) > tolerance.tolerance * baseSize) {
            return std::forward_as_tuple(false, "Block size mismatch.");
        }
    }
    return std::forward_as_tuple(true, "");
}

std::tuple<bool, std::string> ValidationService::ifNodeSizesEqual(const std::vector<Subdomain>& subdomainAssets) {
    if (subdomainAssets.empty()) return std::forward_as_tuple(true, "");
    const size_t baseSize = subdomainAssets[0].nodes.size();
    for (const auto& sub : subdomainAssets) {
        if (sub.nodes.size() != baseSize) return std::forward_as_tuple(false, "Node size mismatch.");
    }
    return std::forward_as_tuple(true, "");
}

std::tuple<bool, std::string> ValidationService::ifNodeOrdersEqual(const std::vector<Subdomain>& subdomainAssets) {
    if (subdomainAssets.empty()) return std::forward_as_tuple(true, "");
    const auto& ref = subdomainAssets[0].nodes;
    for (const auto& sub : subdomainAssets) {
        for (size_t i = 0; i < ref.size(); ++i) {
            if (ref[i]->id != sub.nodes[i]->id) return std::forward_as_tuple(false, "Node order mismatch.");
        }
    }
    return std::forward_as_tuple(true, "");
}

std::tuple<bool, std::string> ValidationService::ifNodePositionsEqual(const std::vector<Subdomain>& subdomainAssets) {
    if (subdomainAssets.empty()) return std::forward_as_tuple(true, "");
    const auto& ref = subdomainAssets[0].nodes;
    for (const auto& sub : subdomainAssets) {
        for (size_t i = 0; i < ref.size(); ++i) {
            const auto& a = ref[i]->pos;
            const auto& b = sub.nodes[i]->pos;
            if (a.x != b.x || a.y != b.y || a.z != b.z) return std::forward_as_tuple(false, "Node position mismatch.");
        }
    }
   return std::forward_as_tuple(true, "");
}

std::tuple<bool, std::string> ValidationService::ifBlocksValid(const Domain& domain, const Tolerance tolerance) {
    const auto& blocks = domain.blocks;

    bool passed;
    std::string errorMessage;

    std::tie(passed, errorMessage) = ifBlockSizesEqual(blocks, tolerance);
    if (!passed) return std::forward_as_tuple(false, errorMessage);

    std::tie(passed, errorMessage) = ifNodeSizesEqual(blocks);
    if (!passed) return std::forward_as_tuple(false, errorMessage);

    std::tie(passed, errorMessage) = ifNodeOrdersEqual(blocks);
    if (!passed) return std::forward_as_tuple(false, errorMessage);

    std::tie(passed, errorMessage) = ifNodePositionsEqual(blocks);
    if (!passed) return std::forward_as_tuple(false, errorMessage);

    return std::forward_as_tuple(true, "");
}

} // namespace partitioner


Args Args::parse(int argc, char* argv[], size_t requiredPositionals) {
    Args args;

    for (int i = 1; i < argc; ++i) {
        std::string token = argv[i];

        if (!token.empty() && token[0] == '-') {
            // violated if option is not included
            if (token.size() == 1) {
                throw std::invalid_argument("Invalid option: '-' is not a valid key.");
            }

            // violated if no next token
            if (i + 1 >= argc) {
                throw std::invalid_argument("Option '" + token + "' requires a value.");
            }

            // violated if next token is an option
            std::string next = argv[i + 1];
            if (!next.empty() && next[0] == '-') {
                throw std::invalid_argument("Option '" + token + "' must have a value.");
            }

            args.options[token] = next;
            ++i;
        } else {
            args.positionals.push_back(token);
        }
    }

    if (args.positionals.size() < requiredPositionals) {
        throw std::invalid_argument(
            "Expected at least " + std::to_string(requiredPositionals) +
            " positional arguments, but got " + std::to_string(args.positionals.size()) + '.');
    }

    return args;
}

std::string Args::getOpt(const std::string& key) const {
    if (options.find(key) == options.end()) {
        throw std::invalid_argument("Option '" + key + "' not found.");
    }
    return options.at(key);
};