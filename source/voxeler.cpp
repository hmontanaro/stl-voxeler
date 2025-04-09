#include "voxeler.h"
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

// Save voxel size and grid dimensions to a separate metadata file
void saveVoxelMetadata(const std::string& baseFilename, double voxelSize, int Nx, int Ny, int Nz, 
    const Vector3& minB, const Vector3& maxB, int offsetX, int offsetY, int offsetZ) {
    // Create metadata filename
    std::string metaFilename = baseFilename + "-metadata.txt";
    std::cout << "Writing voxelization metadata to " << metaFilename << std::endl;

    // Open the file
    std::ofstream metaFile(metaFilename);
    if (!metaFile) {
    std::cerr << "Error: Could not open metadata file for writing." << std::endl;
    return;
    }

    // Write voxel dimensions and other useful information
    metaFile << "voxel_size: " << voxelSize << std::endl;
    metaFile << "grid_dimensions: " << Nx << " " << Ny << " " << Nz << std::endl;
    metaFile << "total_voxels: " << (Nx * Ny * Nz) << std::endl;
    metaFile << "min_bounds: " << minB.x << " " << minB.y << " " << minB.z << std::endl;
    metaFile << "max_bounds: " << maxB.x << " " << maxB.y << " " << maxB.z << std::endl;
    metaFile << "voxel_offset: " << offsetX << " " << offsetY << " " << offsetZ << std::endl;

    metaFile.close();
    std::cout << "Metadata saved successfully." << std::endl;
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
