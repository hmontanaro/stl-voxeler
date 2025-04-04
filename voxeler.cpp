#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <queue>
#include <chrono>
#include <iomanip>

// Vector3 structure for 3D vector representation
struct Vector3 {
    double x, y, z;
    // Default constructor initializes vector to origin (0,0,0)
    Vector3() : x(0), y(0), z(0) {}
    // Constructor with specific coordinates
    Vector3(double a, double b, double c) : x(a), y(b), z(c) {}
    // Vector addition operator
    Vector3 operator+(const Vector3& v) const { return Vector3(x + v.x, y + v.y, z + v.z); }
    // Vector subtraction operator
    Vector3 operator-(const Vector3& v) const { return Vector3(x - v.x, y - v.y, z - v.z); }
    // Scalar multiplication operator
    Vector3 operator*(double s) const { return Vector3(x * s, y * s, z * s); }
};

// Calculate dot product between two vectors (scalar result)
double dot(const Vector3& a, const Vector3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Calculate cross product between two vectors (vector result)
// Returns a vector perpendicular to both input vectors
Vector3 cross(const Vector3& a, const Vector3& b) {
    return Vector3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

// Triangle structure for 3D geometry representation
// Contains three vertices defining the triangle
struct Triangle {
    Vector3 v0, v1, v2;
};

/**
 * Loads triangular mesh from an STL file
 * 
 * @param filename Path to the STL file
 * @param triangles Vector to store the loaded triangles
 * @return true if loading was successful, false otherwise
 */
bool loadSTL(const std::string& filename, std::vector<Triangle>& triangles) {
    // Open file in binary mode first to check format
    std::ifstream file(filename.c_str(), std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file: " << filename << "\n";
        return false;
    }
    
    // Read first 5 bytes to check if it's ASCII (solid) or binary
    // STL ASCII format starts with "solid"
    char header[6] = {0};
    file.read(header, 5);
    header[5] = '\0';  // Null-terminate for string comparison
    file.close();
    
    std::string headerStr(header);
    bool isAscii = (headerStr == "solid");
    
    // If it looks like ASCII, double-check it's not a binary file with "solid" in header
    if (isAscii) {
        file.open(filename.c_str(), std::ios::binary);
        file.seekg(0, std::ios::end);
        size_t fileSize = file.tellg();
        file.close();
        
        // Binary STL with n triangles: 84 + n*50 bytes
        // 84 = 80 byte header + 4 byte triangle count
        // 50 = 12 byte normal + 36 byte vertices + 2 byte attribute
        // If file size matches this pattern, it's likely binary despite "solid" header
        if ((fileSize - 84) % 50 == 0) {
            isAscii = false;
            // TODO: Consider more robust format detection methods
        }
    }
    
    if (isAscii) {
        // ASCII STL format
        file.open(filename.c_str());
        std::string line;
        Triangle tri;
        int vertexCount = 0;
        
        // TODO: This parsing is minimal and might break with non-standard STL files
        // Process the file line by line
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string word;
            iss >> word;
            // When we find a "vertex" line, read the coordinates
            if (word == "vertex") {
                double x, y, z;
                // TODO: Add error checking for incomplete vertex data
                iss >> x >> y >> z;
                if (vertexCount == 0) tri.v0 = Vector3(x, y, z);
                else if (vertexCount == 1) tri.v1 = Vector3(x, y, z);
                else if (vertexCount == 2) tri.v2 = Vector3(x, y, z);
                vertexCount++;
                if (vertexCount == 3) {
                    triangles.push_back(tri);
                    vertexCount = 0;
                }
            }
        }
    } else {
        // Binary STL format
        file.open(filename.c_str(), std::ios::binary);
        
        // Skip 80-byte header
        file.seekg(80);
        
        // Read number of triangles (4-byte unsigned int)
        uint32_t numTriangles;
        file.read(reinterpret_cast<char*>(&numTriangles), sizeof(numTriangles));
        
        // Pre-allocate triangle vector
        triangles.reserve(numTriangles);
        
        // Read each triangle: normal (3 floats), vertices (9 floats), attribute (2 bytes)
        for (uint32_t i = 0; i < numTriangles; i++) {
            float normal[3];
            float vertices[9];
            uint16_t attribute;
            
            // Skip normal vector (12 bytes)
            file.seekg(12, std::ios::cur);
            
            // Read vertices
            file.read(reinterpret_cast<char*>(vertices), sizeof(vertices));
            
            // Skip attribute
            file.seekg(2, std::ios::cur);
            
            // Create triangle
            Triangle tri;
            tri.v0 = Vector3(vertices[0], vertices[1], vertices[2]);
            tri.v1 = Vector3(vertices[3], vertices[4], vertices[5]);
            tri.v2 = Vector3(vertices[6], vertices[7], vertices[8]);
            triangles.push_back(tri);
        }
    }
    
    std::cout << "Loaded STL file in " << (isAscii ? "ASCII" : "binary") << " format\n";
    return !triangles.empty();
}

/**
 * Determines if a triangle from the STL surface mesh intersects with an axis-aligned box (a voxel) using the Separating Axis Theorem (SAT)
 * 
 * This function is critical for the STL to voxel conversion process, as it determines
 * which voxels in the grid are intersected by the triangular mesh.
 * 
 * @param boxcenter Center coordinates of the box (voxel)
 * @param boxhalfsize Half-dimensions of the box (half the voxel size in each direction)
 * @param tri The triangle to test against the box
 * @return true if the triangle and box intersect, false otherwise
 */
bool triBoxOverlap(const Vector3& boxcenter, const Vector3& boxhalfsize, const Triangle& tri) {
    // Translate triangle vertices to be relative to the box center
    Vector3 v0 = tri.v0 - boxcenter;
    Vector3 v1 = tri.v1 - boxcenter;
    Vector3 v2 = tri.v2 - boxcenter;

    // Compute triangle edges
    Vector3 e0 = v1 - v0;  // Edge from v0 to v1
    Vector3 e1 = v2 - v1;  // Edge from v1 to v2
    Vector3 e2 = v0 - v2;  // Edge from v2 to v0

    double min, max, p0, p1, p2, rad;

    // The following tests check for a separating axis between the box and triangle
    // If any separating axis is found, then the shapes don't intersect
    
    // Test nine potential separating axes formed by cross products of
    // triangle edges with the box's primary axes (X, Y, Z)
    
    // Test axis L = e0 × (1,0,0) = (0, e0.z, -e0.y)
    p0 = v0.z * e0.y - v0.y * e0.z;
    p1 = v1.z * e0.y - v1.y * e0.z;
    p2 = v2.z * e0.y - v2.y * e0.z;
    min = std::min(std::min(p0, p1), p2);
    max = std::max(std::max(p0, p1), p2);
    rad = boxhalfsize.y * fabs(e0.z) + boxhalfsize.z * fabs(e0.y);
    if (min > rad || max < -rad) return false;  // Separating axis found

    // Test axis L = e0 × (0,1,0) = (-e0.z, 0, e0.x)
    p0 = -(v0.z * e0.x - v0.x * e0.z);
    p1 = -(v1.z * e0.x - v1.x * e0.z);
    p2 = -(v2.z * e0.x - v2.x * e0.z);
    min = std::min(std::min(p0, p1), p2);
    max = std::max(std::max(p0, p1), p2);
    rad = boxhalfsize.x * fabs(e0.z) + boxhalfsize.z * fabs(e0.x);
    if (min > rad || max < -rad) return false;  // Separating axis found

    // Test axis L = e0 × (0,0,1) = (e0.y, -e0.x, 0)
    p0 = v0.y * e0.x - v0.x * e0.y;
    p1 = v1.y * e0.x - v1.x * e0.y;
    p2 = v2.y * e0.x - v2.x * e0.y;
    min = std::min(std::min(p0, p1), p2);
    max = std::max(std::max(p0, p1), p2);
    rad = boxhalfsize.x * fabs(e0.y) + boxhalfsize.y * fabs(e0.x);
    if (min > rad || max < -rad) return false;  // Separating axis found

    // Test axis L = e1 × (1,0,0) = (0, e1.z, -e1.y)
    p0 = v0.z * e1.y - v0.y * e1.z;
    p1 = v1.z * e1.y - v1.y * e1.z;
    p2 = v2.z * e1.y - v2.y * e1.z;
    min = std::min(std::min(p0, p1), p2);
    max = std::max(std::max(p0, p1), p2);
    rad = boxhalfsize.y * fabs(e1.z) + boxhalfsize.z * fabs(e1.y);
    if (min > rad || max < -rad) return false;  // Separating axis found

    // Test axis L = e1 × (0,1,0) = (-e1.z, 0, e1.x)
    p0 = -(v0.z * e1.x - v0.x * e1.z);
    p1 = -(v1.z * e1.x - v1.x * e1.z);
    p2 = -(v2.z * e1.x - v2.x * e1.z);
    min = std::min(std::min(p0, p1), p2);
    max = std::max(std::max(p0, p1), p2);
    rad = boxhalfsize.x * fabs(e1.z) + boxhalfsize.z * fabs(e1.x);
    if (min > rad || max < -rad) return false;  // Separating axis found

    // Test axis L = e1 × (0,0,1) = (e1.y, -e1.x, 0)
    p0 = v0.y * e1.x - v0.x * e1.y;
    p1 = v1.y * e1.x - v1.x * e1.y;
    p2 = v2.y * e1.x - v2.x * e1.y;
    min = std::min(std::min(p0, p1), p2);
    max = std::max(std::max(p0, p1), p2);
    rad = boxhalfsize.x * fabs(e1.y) + boxhalfsize.y * fabs(e1.x);
    if (min > rad || max < -rad) return false;  // Separating axis found

    // Test axis L = e2 × (1,0,0) = (0, e2.z, -e2.y)
    p0 = v0.z * e2.y - v0.y * e2.z;
    p1 = v1.z * e2.y - v1.y * e2.z;
    p2 = v2.z * e2.y - v2.y * e2.z;
    min = std::min(std::min(p0, p1), p2);
    max = std::max(std::max(p0, p1), p2);
    rad = boxhalfsize.y * fabs(e2.z) + boxhalfsize.z * fabs(e2.y);
    if (min > rad || max < -rad) return false;  // Separating axis found

    // Test axis L = e2 × (0,1,0) = (-e2.z, 0, e2.x)
    p0 = -(v0.z * e2.x - v0.x * e2.z);
    p1 = -(v1.z * e2.x - v1.x * e2.z);
    p2 = -(v2.z * e2.x - v2.x * e2.z);
    min = std::min(std::min(p0, p1), p2);
    max = std::max(std::max(p0, p1), p2);
    rad = boxhalfsize.x * fabs(e2.z) + boxhalfsize.z * fabs(e2.x);
    if (min > rad || max < -rad) return false;  // Separating axis found

    // Test axis L = e2 × (0,0,1) = (e2.y, -e2.x, 0)
    p0 = v0.y * e2.x - v0.x * e2.y;
    p1 = v1.y * e2.x - v1.x * e2.y;
    p2 = v2.y * e2.x - v2.x * e2.y;
    min = std::min(std::min(p0, p1), p2);
    max = std::max(std::max(p0, p1), p2);
    rad = boxhalfsize.x * fabs(e2.y) + boxhalfsize.y * fabs(e2.x);
    if (min > rad || max < -rad) return false;  // Separating axis found

    // Test the three box axes (X, Y, Z) - box face normals
    double minTri, maxTri;
    
    // Test X-axis separation
    minTri = std::min(std::min(v0.x, v1.x), v2.x);
    maxTri = std::max(std::max(v0.x, v1.x), v2.x);
    if (minTri > boxhalfsize.x || maxTri < -boxhalfsize.x) return false;
    
    // Test Y-axis separation
    minTri = std::min(std::min(v0.y, v1.y), v2.y);
    maxTri = std::max(std::max(v0.y, v1.y), v2.y);
    if (minTri > boxhalfsize.y || maxTri < -boxhalfsize.y) return false;
    
    // Test Z-axis separation
    minTri = std::min(std::min(v0.z, v1.z), v2.z);
    maxTri = std::max(std::max(v0.z, v1.z), v2.z);
    if (minTri > boxhalfsize.z || maxTri < -boxhalfsize.z) return false;

    // Test separation by the triangle face normal
    Vector3 normal = cross(e0, e1);  // Triangle normal
    double d = -dot(normal, v0);     // Distance from origin to triangle plane
    
    // Find box vertices with min/max projection along the normal
    Vector3 vmin, vmax;
    vmin.x = (normal.x > 0.0) ? -boxhalfsize.x : boxhalfsize.x;
    vmax.x = (normal.x > 0.0) ? boxhalfsize.x : -boxhalfsize.x;
    vmin.y = (normal.y > 0.0) ? -boxhalfsize.y : boxhalfsize.y;
    vmax.y = (normal.y > 0.0) ? boxhalfsize.y : -boxhalfsize.y;
    vmin.z = (normal.z > 0.0) ? -boxhalfsize.z : boxhalfsize.z;
    vmax.z = (normal.z > 0.0) ? boxhalfsize.z : -boxhalfsize.z;
    
    // Test if the box is completely on the negative side of the triangle plane
    if (dot(normal, vmin) + d > 0.0) return false;
    
    // Test if the box intersects the triangle plane
    if (dot(normal, vmax) + d >= 0.0) return true;
    
    // No separating axis found, the triangle and box must intersect
    return false;
}

struct Coord {
    int i, j, k;
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: voxelizer <input.stl> [voxel_size] [fill]\n";
        return 1;
    }
    std::string filename = argv[1];
    double voxelSize = 1.0;
    if (argc >= 3) {
        voxelSize = atof(argv[2]);
        if (voxelSize <= 0) voxelSize = 1.0;
    }

    bool fillInterior = false;
    if (argc >= 4) {
        std::string flag = argv[3];
        if (flag == "fill" || flag == "--fill")
            fillInterior = true;
    }

    std::vector<Triangle> triangles;
    if (!loadSTL(filename, triangles))
        return 1;
    std::cout << "Loaded " << triangles.size() << " triangles.\n";

    Vector3 minB(std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max());
    Vector3 maxB(-std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max());
    for (const auto& tri : triangles) {
        for (const auto& v : { tri.v0, tri.v1, tri.v2 }) {
            if (v.x < minB.x) minB.x = v.x;
            if (v.y < minB.y) minB.y = v.y;
            if (v.z < minB.z) minB.z = v.z;
            if (v.x > maxB.x) maxB.x = v.x;
            if (v.y > maxB.y) maxB.y = v.y;
            if (v.z > maxB.z) maxB.z = v.z;
        }
    }
    double padding = voxelSize;
    minB = minB - Vector3(padding, padding, padding);
    maxB = maxB + Vector3(padding, padding, padding);

    int Nx = static_cast<int>(std::ceil((maxB.x - minB.x) / voxelSize));
    int Ny = static_cast<int>(std::ceil((maxB.y - minB.y) / voxelSize));
    int Nz = static_cast<int>(std::ceil((maxB.z - minB.z) / voxelSize));
    int Nxyz = Nx * Ny * Nz;
    std::cout << "Grid dimensions: " << Nx << " x " << Ny << " x " << Nz << "\n";
    std::cout << "Grid size: " << Nxyz << "\n";

    std::vector<bool> voxels(Nxyz, false);
    
    // Calculate and print memory usage of the voxel array
    size_t voxelCount = voxels.size();
    size_t memoryUsageBytes = voxelCount * sizeof(bool);
    double memoryUsageMB = memoryUsageBytes / (1024.0 * 1024.0);
    
    std::cout << "Voxel array memory usage: " 
              << memoryUsageBytes << " bytes (" 
              << std::fixed << std::setprecision(2) << memoryUsageMB << " MB)" 
              << std::endl;

    auto index = [=](int i, int j, int k) -> int {
        return i + Nx * (j + Ny * k);
    };

    std::cout << "Voxeling triangles\n";
    int triangleCount = triangles.size();
    int processedCount = 0;
    int lastPercentage = -1;
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Counter for occupied voxels
    int occupiedVoxelCount = 0;
    
    for (const auto& tri : triangles) {
        Vector3 triMin(std::min({ tri.v0.x, tri.v1.x, tri.v2.x }),
            std::min({ tri.v0.y, tri.v1.y, tri.v2.y }),
            std::min({ tri.v0.z, tri.v1.z, tri.v2.z }));
        Vector3 triMax(std::max({ tri.v0.x, tri.v1.x, tri.v2.x }),
            std::max({ tri.v0.y, tri.v1.y, tri.v2.y }),
            std::max({ tri.v0.z, tri.v1.z, tri.v2.z }));
        int i0 = std::max(0, static_cast<int>(std::floor((triMin.x - minB.x) / voxelSize)));
        int j0 = std::max(0, static_cast<int>(std::floor((triMin.y - minB.y) / voxelSize)));
        int k0 = std::max(0, static_cast<int>(std::floor((triMin.z - minB.z) / voxelSize)));
        int i1 = std::min(Nx - 1, static_cast<int>(std::floor((triMax.x - minB.x) / voxelSize)));
        int j1 = std::min(Ny - 1, static_cast<int>(std::floor((triMax.y - minB.y) / voxelSize)));
        int k1 = std::min(Nz - 1, static_cast<int>(std::floor((triMax.z - minB.z) / voxelSize)));

        for (int k = k0; k <= k1; ++k) {
            for (int j = j0; j <= j1; ++j) {
                for (int i = i0; i <= i1; ++i) {
                    Vector3 voxelCenter(minB.x + (i + 0.5) * voxelSize,
                        minB.y + (j + 0.5) * voxelSize,
                        minB.z + (k + 0.5) * voxelSize);
                    Vector3 halfSize(voxelSize / 2.0, voxelSize / 2.0, voxelSize / 2.0);
                    if (triBoxOverlap(voxelCenter, halfSize, tri)) {
                        if (!voxels[index(i, j, k)]) {
                            // Only count if this voxel wasn't already marked
                            occupiedVoxelCount++;
                        }
                        voxels[index(i, j, k)] = true;
                    }
                }
            }
        }
        
        // Update progress every 10%
        processedCount++;
        int currentPercentage = (processedCount * 100) / triangleCount;
        if (currentPercentage / 10 > lastPercentage / 10) {
            std::cout << " " << (currentPercentage / 10) * 10 << "%" << std::endl;
            lastPercentage = currentPercentage;
        }
    }
    
    std::cout << "Occupied: " << occupiedVoxelCount << "/" << Nxyz << " (" << (occupiedVoxelCount * 100.0 / Nxyz) << "%)" << std::endl;
    
    // Ensure 100% is printed at the end if not already printed
    if (lastPercentage < 100) {
        std::cout << "Progress: 100%" << std::endl;
    }

	auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "  Completed in " << duration.count() / 1000.0 << " seconds" << std::endl;

    if (fillInterior) {
        std::cout << "Filling interior volume\n";

        // Ray casting fill: process each slice along z and each row along y.
        int totalSlices = Nz;
        int processedSlices = 0;
        int lastReportedPercentage = 0;
        int interiorVoxelsFilled = 0;
        
        auto startTime = std::chrono::high_resolution_clock::now();
        
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                bool inside = false;
                for (int i = 0; i < Nx; ++i) {
                    int idx = index(i, j, k);
                    // Toggle "inside" when entering a new boundary segment.
                    if (voxels[idx]) {
                        if (i == 0 || !voxels[index(i - 1, j, k)]) {
                            inside = !inside;
                        }
                    }
                    else {
                        if (inside) {
                            voxels[idx] = true;
                            interiorVoxelsFilled++;
                        }
                    }
                }
            }
            
            // Update progress
            processedSlices++;
            int currentPercentage = (processedSlices * 100) / totalSlices;
            if (currentPercentage / 10 > lastReportedPercentage / 10) {
                lastReportedPercentage = currentPercentage;
                std::cout << " " << currentPercentage << "%" << std::endl;
            }
        }
        
        std::cout << "Interior voxels filled: " << interiorVoxelsFilled << std::endl;
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        std::cout << "  Completed in " << duration.count() / 1000.0 << " seconds" << std::endl;

            std::cout << "  Voxel array memory usage: " 
              << memoryUsageBytes << " bytes (" 
              << std::fixed << std::setprecision(2) << memoryUsageMB << " MB)" 
              << std::endl;
    }
	
    std::cout << "Counting occupied voxels..." << std::endl;
    auto countStartTime = std::chrono::high_resolution_clock::now();
    
    int count = 0;
    for (bool v : voxels)
        if (v) count++;
    
    auto countEndTime = std::chrono::high_resolution_clock::now();
    auto countDuration = std::chrono::duration_cast<std::chrono::milliseconds>(countEndTime - countStartTime);
    std::cout << "Occupied voxels: " << count << " (counted in " 
              << countDuration.count() / 1000.0 << " seconds)" << std::endl;

    // Create output filename based on input filename without extension
    std::string baseFilename = filename;
    size_t lastDot = baseFilename.find_last_of(".");
    if (lastDot != std::string::npos) {
        baseFilename = baseFilename.substr(0, lastDot);
    }
    std::string outputFilename = baseFilename + "-voxeled.txt";
    
    // std::ofstream out(outputFilename);
    // for (int k = 0; k < Nz; ++k)
    //     for (int j = 0; j < Ny; ++j)
    //         for (int i = 0; i < Nx; ++i)
    //             if (voxels[index(i, j, k)]) {
    //                 Vector3 center(minB.x + (i + 0.5) * voxelSize,
    //                     minB.y + (j + 0.5) * voxelSize,
    //                     minB.z + (k + 0.5) * voxelSize);
    //                 out << center.x << " " << center.y << " " << center.z << "\n";
    //             }
    // out.close();

    // Text format with optimized buffering
    std::ofstream out(outputFilename);
    
    // Set larger buffer size (8MB)
    const size_t BUFFER_SIZE = 8 * 1024 * 1024;
    char* buffer = new char[BUFFER_SIZE];
    out.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);
    
    // Disable synchronization with C standard streams
    out.sync_with_stdio(false);
    
    // Use a string stream for efficient string building
    std::stringstream ss;
    ss.sync_with_stdio(false);
    
    // Process in chunks to avoid excessive memory usage
    const int CHUNK_SIZE = 1000000; // Process 1M voxels at a time
    std::vector<Vector3> centers;
    centers.reserve(CHUNK_SIZE);
    
    // First pass: collect voxel centers
    auto writeStartTime = std::chrono::high_resolution_clock::now();
    int totalVoxels = Nx * Ny * Nz;
    int processedVoxels = 0;
    lastPercentage = -1;
    
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                // Update progress every 1%
                ++processedVoxels;
                int currentPercentage = (processedVoxels * 100) / totalVoxels;
                if (currentPercentage > lastPercentage) {
                    std::cout << "\rWriting voxels: " << currentPercentage << "% complete" << std::flush;
                    lastPercentage = currentPercentage;
                }
                
                if (voxels[index(i, j, k)]) {
                    Vector3 center(minB.x + (i + 0.5) * voxelSize,
                                    minB.y + (j + 0.5) * voxelSize,
                                    minB.z + (k + 0.5) * voxelSize);
                    centers.push_back(center);
                    
                    // Write chunk if buffer is full
                    if (centers.size() >= CHUNK_SIZE) {
                        for (const auto& c : centers) {
                            ss << c.x << " " << c.y << " " << c.z << "\n";
                        }
                        out << ss.str();
                        ss.str(""); // Clear the string stream
                        centers.clear();
                    }
                }
            }
        }
    }
    
    // Write any remaining centers
    if (!centers.empty()) {
        for (const auto& c : centers) {
            ss << c.x << " " << c.y << " " << c.z << "\n";
        }
        out << ss.str();
    }
    
    out.close();
    delete[] buffer;
    
    auto writeEndTime = std::chrono::high_resolution_clock::now();
    auto writeDuration = std::chrono::duration_cast<std::chrono::milliseconds>(writeEndTime - writeStartTime);
    std::cout << "\nVoxel centers saved to " << outputFilename << " in " 
              << writeDuration.count() / 1000.0 << " seconds\n";

    return 0;
}
