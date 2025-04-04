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

/**
 * Coordinate structure for voxel grid indexing
 * Stores integer coordinates in the voxel grid
 */
struct Coord {
    int i, j, k;  // i, j, k represent indices along the x, y, z axes respectively
};

/**
 * Main function for STL voxelization process
 * Handles command-line arguments, loads STL file, and converts it to voxels
 */
int main(int argc, char* argv[]) {
    // Validate command-line arguments
    if (argc < 2) {
        std::cerr << "Usage: voxelizer <input.stl> [voxel_size] [fill]\n";
        return 1;
    }
    
    // Extract input filename from command-line arguments
    std::string filename = argv[1];
    
    // Parse optional voxel size parameter (default: 1.0)
    // Determines the resolution of the voxelization - smaller values create more detailed output
    double voxelSize = 1.0;
    if (argc >= 3) {
        voxelSize = atof(argv[2]);
        // Safety check to avoid invalid voxel sizes
        if (voxelSize <= 0) voxelSize = 1.0;
    }

    // Parse optional fill flag which determines whether to fill the interior volume
    // If enabled, the algorithm will mark voxels inside the mesh as occupied
    bool fillInterior = false;
    if (argc >= 4) {
        std::string flag = argv[3];
        if (flag == "fill" || flag == "--fill")
            fillInterior = true;
    }

    // Load STL mesh data from the input file
    // Returns error if file loading fails or format is invalid
    std::vector<Triangle> triangles;
    if (!loadSTL(filename, triangles))
        return 1;
    std::cout << "Loaded " << triangles.size() << " triangles.\n";

    // Calculate the bounding box of the mesh to determine voxel grid dimensions
    // First set to maximum/minimum possible values, then update with actual bounds
    Vector3 minB(std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max());
    Vector3 maxB(-std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max());
    
    // Iterate through all triangles and their vertices to find min/max coordinates
    for (const auto& tri : triangles) {
        for (const auto& v : { tri.v0, tri.v1, tri.v2 }) {
            // Update minimum bounds
            if (v.x < minB.x) minB.x = v.x;
            if (v.y < minB.y) minB.y = v.y;
            if (v.z < minB.z) minB.z = v.z;
            // Update maximum bounds
            if (v.x > maxB.x) maxB.x = v.x;
            if (v.y > maxB.y) maxB.y = v.y;
            if (v.z > maxB.z) maxB.z = v.z;
        }
    }
    
    // Add padding to ensure the entire mesh is inside the voxel grid
    // This prevents clipping of the mesh at the boundaries
    double padding = voxelSize;
    minB = minB - Vector3(padding, padding, padding);
    maxB = maxB + Vector3(padding, padding, padding);

    // Calculate voxel grid dimensions based on the bounding box and voxel size
    // ceil() ensures we have enough voxels to cover the entire bounding box
    int Nx = static_cast<int>(std::ceil((maxB.x - minB.x) / voxelSize));
    int Ny = static_cast<int>(std::ceil((maxB.y - minB.y) / voxelSize));
    int Nz = static_cast<int>(std::ceil((maxB.z - minB.z) / voxelSize));
    int Nxyz = Nx * Ny * Nz; // Total number of voxels in the grid
    std::cout << "Grid dimensions: " << Nx << " x " << Ny << " x " << Nz << "\n";
    std::cout << "Grid size: " << Nxyz << "\n";

    // Initialize voxel grid as a flat 1D array for better memory efficiency
    // All voxels are initially empty (false)
    std::vector<bool> voxels(Nxyz, false);
    
    // Calculate and display memory usage of the voxel array
    // This helps users understand the resource requirements of the voxelization
    size_t voxelCount = voxels.size();
    size_t memoryUsageBytes = voxelCount * sizeof(bool);
    double memoryUsageMB = memoryUsageBytes / (1024.0 * 1024.0);
    
    std::cout << "Voxel array memory usage: " 
              << memoryUsageBytes << " bytes (" 
              << std::fixed << std::setprecision(2) << memoryUsageMB << " MB)" 
              << std::endl;

    // Define a lambda function to convert 3D coordinates to 1D array index
    // This enables efficient access to the voxel grid using i,j,k coordinates
    auto index = [=](int i, int j, int k) -> int {
        return i + Nx * (j + Ny * k);
    };

    std::cout << "Voxeling triangles\n";
    int triangleCount = triangles.size();
    int processedCount = 0;
    int lastPercentage = -1;
    
    // Record start time for performance measurement
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Track how many voxels are marked as occupied
    int occupiedVoxelCount = 0;
    
    // Process each triangle to identify intersecting voxels
    for (const auto& tri : triangles) {
        // Calculate bounding box of the triangle to limit voxel checks
        // This optimization avoids checking every voxel for every triangle
        Vector3 triMin(std::min({ tri.v0.x, tri.v1.x, tri.v2.x }),
            std::min({ tri.v0.y, tri.v1.y, tri.v2.y }),
            std::min({ tri.v0.z, tri.v1.z, tri.v2.z }));
        Vector3 triMax(std::max({ tri.v0.x, tri.v1.x, tri.v2.x }),
            std::max({ tri.v0.y, tri.v1.y, tri.v2.y }),
            std::max({ tri.v0.z, tri.v1.z, tri.v2.z }));
        
        // Convert triangle bounds to voxel grid indices
        // These determine the subset of voxels to check for intersection
        int i0 = std::max(0, static_cast<int>(std::floor((triMin.x - minB.x) / voxelSize)));
        int j0 = std::max(0, static_cast<int>(std::floor((triMin.y - minB.y) / voxelSize)));
        int k0 = std::max(0, static_cast<int>(std::floor((triMin.z - minB.z) / voxelSize)));
        int i1 = std::min(Nx - 1, static_cast<int>(std::floor((triMax.x - minB.x) / voxelSize)));
        int j1 = std::min(Ny - 1, static_cast<int>(std::floor((triMax.y - minB.y) / voxelSize)));
        int k1 = std::min(Nz - 1, static_cast<int>(std::floor((triMax.z - minB.z) / voxelSize)));

        // Check each voxel within the triangle's bounding box for intersection
        for (int k = k0; k <= k1; ++k) {
            for (int j = j0; j <= j1; ++j) {
                for (int i = i0; i <= i1; ++i) {
                    // Calculate center of the current voxel in world coordinates
                    Vector3 voxelCenter(minB.x + (i + 0.5) * voxelSize,
                        minB.y + (j + 0.5) * voxelSize,
                        minB.z + (k + 0.5) * voxelSize);
                    // Define half-size of the voxel cube
                    Vector3 halfSize(voxelSize / 2.0, voxelSize / 2.0, voxelSize / 2.0);
                    
                    // Test if the triangle intersects this voxel using SAT algorithm
                    if (triBoxOverlap(voxelCenter, halfSize, tri)) {
                        if (!voxels[index(i, j, k)]) {
                            // Only increment if this voxel wasn't already marked
                            // This prevents double-counting when multiple triangles intersect a voxel
                            occupiedVoxelCount++;
                        }
                        voxels[index(i, j, k)] = true;
                    }
                }
            }
        }
        
        // Update progress display every 10%
        // This helps users monitor the voxelization process for large models
        processedCount++;
        int currentPercentage = (processedCount * 100) / triangleCount;
        if (currentPercentage / 10 > lastPercentage / 10) {
            std::cout << " " << (currentPercentage / 10) * 10 << "%" << std::endl;
            lastPercentage = currentPercentage;
        }
    }
    
    std::cout << "Occupied: " << occupiedVoxelCount << "/" << Nxyz << " (" << (occupiedVoxelCount * 100.0 / Nxyz) << "%)" << std::endl;
    
    // Ensure 100% progress is reported if not already printed
    if (lastPercentage < 100) {
        std::cout << "Progress: 100%" << std::endl;
    }

    // Output the total processing time for the voxelization step
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "  Completed in " << duration.count() / 1000.0 << " seconds" << std::endl;
    // Interior filling using flood fill algorithm
    // This marks all voxels inside the mesh boundary as occupied
    if (fillInterior) {
        std::cout << "Filling interior using flood fill algorithm\n";
        
        // Start timing
        auto startTime = std::chrono::high_resolution_clock::now();
        
        // Create visited grid to track the flood fill
        std::vector<bool> visited(Nxyz, false);
        
        // Initialize the exterior voxel count
        int exteriorVoxels = 0;
        
        // Function to check if a voxel is valid and can be filled
        auto isValid = [&](int i, int j, int k) -> bool {
            // Check if indices are within bounds
            if (i < 0 || i >= Nx || j < 0 || j >= Ny || k < 0 || k >= Nz)
                return false;
                
            int idx = index(i, j, k);
            // Valid if not already occupied by mesh and not already visited
            return !voxels[idx] && !visited[idx];
        };
        
        // Use BFS (Breadth-First Search) for flood fill
        std::queue<Coord> queue;
        
        // Start with boundary voxels (6 faces of the volume)
        // This ensures we begin the fill from outside the mesh
        
        // Add voxels on the X-boundary planes
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i : {0, Nx - 1}) {
                    int idx = index(i, j, k);
                    if (!voxels[idx] && !visited[idx]) {
                        queue.push({i, j, k});
                        visited[idx] = true;
                        exteriorVoxels++;
                    }
                }
            }
        }
        
        // Add voxels on the Y-boundary planes
        for (int k = 0; k < Nz; ++k) {
            for (int i = 1; i < Nx - 1; ++i) {
                for (int j : {0, Ny - 1}) {
                    int idx = index(i, j, k);
                    if (!voxels[idx] && !visited[idx]) {
                        queue.push({i, j, k});
                        visited[idx] = true;
                        exteriorVoxels++;
                    }
                }
            }
        }
        
        // Add voxels on the Z-boundary planes
        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                for (int k : {0, Nz - 1}) {
                    int idx = index(i, j, k);
                    if (!voxels[idx] && !visited[idx]) {
                        queue.push({i, j, k});
                        visited[idx] = true;
                        exteriorVoxels++;
                    }
                }
            }
        }
        
        // Define the 6 directions for flood fill (3D neighbors)
        const int dx[] = {1, -1, 0, 0, 0, 0};
        const int dy[] = {0, 0, 1, -1, 0, 0};
        const int dz[] = {0, 0, 0, 0, 1, -1};
        
        // Process the queue for BFS flood fill
        int lastPercentage = 0;
        size_t totalVoxels = voxels.size();
        size_t processedVoxels = 0;
        
        while (!queue.empty()) {
            Coord current = queue.front();
            queue.pop();
            
            // Progress tracking
            processedVoxels++;
            if (processedVoxels % 100000 == 0 || queue.empty()) {
                int currentPercentage = (processedVoxels * 100) / totalVoxels;
                if (currentPercentage != lastPercentage) {
                    std::cout << "\rFlood fill progress: " << currentPercentage << "% " 
                              << "Queue size: " << queue.size() 
                              << "     " << std::flush;
                    lastPercentage = currentPercentage;
                }
            }
            
            // Try all 6 neighboring directions
            for (int dir = 0; dir < 6; ++dir) {
                int ni = current.i + dx[dir];
                int nj = current.j + dy[dir];
                int nk = current.k + dz[dir];
                
                if (isValid(ni, nj, nk)) {
                    int idx = index(ni, nj, nk);
                    visited[idx] = true;
                    queue.push({ni, nj, nk});
                    exteriorVoxels++;
                }
            }
        }
        
        std::cout << "\nExterior voxels identified: " << exteriorVoxels << std::endl;
        
        // Now identify and mark interior voxels
        int interiorVoxelsFilled = 0;
        
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    int idx = index(i, j, k);
                    // If not visited (exterior) and not already part of mesh surface
                    if (!visited[idx] && !voxels[idx]) {
                        voxels[idx] = true; // Mark as interior
                        interiorVoxelsFilled++;
                    }
                }
            }
        }
        
        std::cout << "Interior voxels filled: " << interiorVoxelsFilled << std::endl;
        
        // Report performance metrics for flood fill
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        std::cout << "  Completed in " << duration.count() / 1000.0 << " seconds" << std::endl;
    }
	
    //// Interior filling using ray-casting algorithm
    //// Marks voxels inside the mesh as occupied
    //if (fillInterior) {
    //    std::cout << "Filling interior volume\n";

    //    // Ray casting fill: process each slice along z and each row along y
    //    // This implements a 2D scan-line fill approach on each z-slice
    //    int totalSlices = Nz;
    //    int processedSlices = 0;
    //    int lastReportedPercentage = 0;
    //    int interiorVoxelsFilled = 0;
    //    
    //    auto startTime = std::chrono::high_resolution_clock::now();
    //    
    //    // Process each z-slice of the voxel grid
    //    for (int k = 0; k < Nz; ++k) {
    //        // Process each row in the slice
    //        for (int j = 0; j < Ny; ++j) {
    //            // Reset parity for each new ray
    //            bool inside = false;
    //            bool wasOccupied = false;
    //            
    //            // Scan along the x-axis (i-coordinate)
    //            for (int i = 0; i < Nx; ++i) {
    //                int idx = index(i, j, k);
    //                bool isOccupied = voxels[idx];
    //                
    //                // Detect transition (boundary crossing)
    //                if (isOccupied != wasOccupied) {
    //                    inside = !inside;
    //                    wasOccupied = isOccupied;
    //                }
    //                
    //                // Fill empty voxels that are inside
    //                if (!isOccupied && inside) {
    //                    voxels[idx] = true;
    //                    interiorVoxelsFilled++;
    //                }
    //            }
    //        }
    //        
    //        // Update progress reporting for interior filling
    //        processedSlices++;
    //        int currentPercentage = (processedSlices * 100) / totalSlices;
    //        if (currentPercentage / 10 > lastReportedPercentage / 10) {
    //            lastReportedPercentage = currentPercentage;
    //            std::cout << " " << currentPercentage << "%" << std::endl;
    //        }
    //    }
    //    
    //    std::cout << "Interior voxels filled: " << interiorVoxelsFilled << std::endl;
    //    
    //    // Report performance metrics for the interior filling operation
    //    auto endTime = std::chrono::high_resolution_clock::now();
    //    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    //    std::cout << "  Completed in " << duration.count() / 1000.0 << " seconds" << std::endl;

    //    // Remind the user about memory usage after interior filling
    //    std::cout << "  Voxel array memory usage: " 
    //          << memoryUsageBytes << " bytes (" 
    //          << std::fixed << std::setprecision(2) << memoryUsageMB << " MB)" 
    //          << std::endl;
    //}
    // Count occupied voxels for final statistics
    std::cout << "Counting occupied voxels..." << std::endl;
    auto countStartTime = std::chrono::high_resolution_clock::now();
    
    int count = 0;
    for (bool v : voxels)
        if (v) count++;
    
    auto countEndTime = std::chrono::high_resolution_clock::now();
    auto countDuration = std::chrono::duration_cast<std::chrono::milliseconds>(countEndTime - countStartTime);
    std::cout << "Occupied voxels: " << count << " (counted in " 
              << countDuration.count() / 1000.0 << " seconds)" << std::endl;

    // Create output filename based on input filename
    // Appends "-voxeled.txt" to the original filename without extension
    std::string baseFilename = filename;
    size_t lastDot = baseFilename.find_last_of(".");
    if (lastDot != std::string::npos) {
        baseFilename = baseFilename.substr(0, lastDot);
    }
    std::string outputFilename = baseFilename + "-voxeled.txt";
    
    // Note: The simple approach below is commented out as it's less efficient
    // for large models due to frequent I/O operations
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

    // Optimized file writing with buffering for better performance
    // This approach significantly reduces I/O overhead for large models
    std::ofstream out(outputFilename);
    
    // Set larger buffer size (8MB) to reduce system calls
    // This improves performance by batching write operations
    const size_t BUFFER_SIZE = 8 * 1024 * 1024;
    char* buffer = new char[BUFFER_SIZE];
    out.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);
    
    // Disable synchronization with C standard streams for better performance
    out.sync_with_stdio(false);
    
    // Use a string stream for efficient string building before writing to file
    std::stringstream ss;
    ss.sync_with_stdio(false);
    
    // Process in chunks to avoid excessive memory usage
    // This allows handling very large models without running out of memory
    const int CHUNK_SIZE = 1000000; // Process 1M voxels at a time
    std::vector<Vector3> centers;
    centers.reserve(CHUNK_SIZE);
    
    // First pass: collect voxel centers in batches and write to file
    auto writeStartTime = std::chrono::high_resolution_clock::now();
    int totalVoxels = Nx * Ny * Nz;
    int processedVoxels = 0;
    lastPercentage = -1;
    
    // Iterate through all voxels in the grid
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
                
                // For occupied voxels, calculate center coordinates and store for output
                if (voxels[index(i, j, k)]) {
                    Vector3 center(minB.x + (i + 0.5) * voxelSize,
                                    minB.y + (j + 0.5) * voxelSize,
                                    minB.z + (k + 0.5) * voxelSize);
                    centers.push_back(center);
                    
                    // Write chunk if buffer is full to manage memory usage
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
    
    // Write any remaining centers at the end
    if (!centers.empty()) {
        for (const auto& c : centers) {
            ss << c.x << " " << c.y << " " << c.z << "\n";
        }
        out << ss.str();
    }
    
    // Clean up resources
    out.close();
    delete[] buffer;
    
    // Report file writing performance
    auto writeEndTime = std::chrono::high_resolution_clock::now();
    auto writeDuration = std::chrono::duration_cast<std::chrono::milliseconds>(writeEndTime - writeStartTime);
    std::cout << "\nVoxel centers saved to " << outputFilename << " in " 
              << writeDuration.count() / 1000.0 << " seconds\n";

    return 0;
}
