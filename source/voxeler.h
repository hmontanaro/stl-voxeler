#ifndef VOXELER_H
#define VOXELER_H

#include <vector>
#include <string>

/**
 * Vector3 structure for 3D vector representation
 */
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

/**
 * Calculate dot product between two vectors (scalar result)
 *
 * @param a First vector
 * @param b Second vector
 * @return Scalar dot product result
 */
double dot(const Vector3& a, const Vector3& b);

/**
 * Calculate cross product between two vectors (vector result)
 * Returns a vector perpendicular to both input vectors
 *
 * @param a First vector
 * @param b Second vector
 * @return Vector cross product result
 */
Vector3 cross(const Vector3& a, const Vector3& b);

/**
 * Triangle structure for 3D geometry representation
 * Contains three vertices defining the triangle
 */
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
bool loadSTL(const std::string& filename, std::vector<Triangle>& triangles);

/**
 * Save voxel size and grid dimensions to a separate metadata file
 *
 * @param baseFilename Base name for the output file
 * @param voxelSize Size of each voxel
 * @param Nx Number of voxels in X direction
 * @param Ny Number of voxels in Y direction
 * @param Nz Number of voxels in Z direction
 * @param minB Minimum bounds of the voxel grid
 * @param maxB Maximum bounds of the voxel grid
 * @param offsetX Optional X offset in voxel units from global origin
 * @param offsetY Optional Y offset in voxel units from global origin
 * @param offsetZ Optional Z offset in voxel units from global origin
 */
void saveVoxelMetadata(const std::string& baseFilename, double voxelSize, int Nx, int Ny, int Nz,
    const Vector3& minB, const Vector3& maxB, int offsetX = 0, int offsetY = 0, int offsetZ = 0);

/**
 * Determines if a triangle from the STL surface mesh intersects with an axis-aligned box (a voxel)
 * using the Separating Axis Theorem (SAT)
 *
 * This function is critical for the STL to voxel conversion process, as it determines
 * which voxels in the grid are intersected by the triangular mesh.
 *
 * @param boxcenter Center coordinates of the box (voxel)
 * @param boxhalfsize Half-dimensions of the box (half the voxel size in each direction)
 * @param tri The triangle to test against the box
 * @return true if the triangle and box intersect, false otherwise
 */
bool triBoxOverlap(const Vector3& boxcenter, const Vector3& boxhalfsize, const Triangle& tri);

/**
 * Voxel index structure for voxel grid indexing
 * Stores integer coordinates in the voxel grid
 */
struct VoxelIdx {
    int i, j, k;  // i, j, k represent indices along the x, y, z axes respectively
};

#endif // VOXELER_H