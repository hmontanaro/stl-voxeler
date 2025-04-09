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

# include "voxeler.h"

/**
 * Main function for STL voxelization process
 * Handles command-line arguments, loads STL file, and converts it to voxels
 */
int main(int argc, char* argv[]) {
    // Validate command-line arguments
    if (argc < 2) {
        std::cerr << "Usage: voxelizer <input.stl> [voxel_size] [options]\n";
        std::cerr << "Options:\n";
        std::cerr << "  --fill      Fill interior volume\n";
        std::cerr << "  --dense     Output in dense format\n";
        std::cerr << "  --align     Align voxels to global grid (0,0,0 origin)\n";
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
    bool denseOutput = false;
    bool gridAligned = false;

    // Process all optional arguments
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "fill" || arg == "--fill") {
            fillInterior = true;
        } else if (arg == "--dense") {
            denseOutput = true;
        } else if (arg == "--align") {
            gridAligned = true;
        }
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
    
    // Print the bounding box information
    std::cout << "\nSTL Bounding Box Information:" << std::endl;
    std::cout << "  X: " << minB.x << " to " << maxB.x << " (width: " << maxB.x - minB.x << ")" << std::endl;
    std::cout << "  Y: " << minB.y << " to " << maxB.y << " (depth: " << maxB.y - minB.y << ")" << std::endl;
    std::cout << "  Z: " << minB.z << " to " << maxB.z << " (height: " << maxB.z - minB.z << ")" << std::endl;
    std::cout << "  Volume: " << (maxB.x - minB.x) * (maxB.y - minB.y) * (maxB.z - minB.z) << std::endl;
    
    // Check for negative coordinates and issue a warning
    if (minB.x < 0 || minB.y < 0 || minB.z < 0) {
        std::cout << "\nWARNING: Model extends into negative coordinate space!" << std::endl;
        if (minB.x < 0) std::cout << "  X minimum is negative: " << minB.x << std::endl;
        if (minB.y < 0) std::cout << "  Y minimum is negative: " << minB.y << std::endl;
        if (minB.z < 0) std::cout << "  Z minimum is negative: " << minB.z << std::endl;
    }
    
    // Calculate voxel offset from global origin
    // This represents how many voxels away from (0,0,0) our grid starts
    int offsetX = 0, offsetY = 0, offsetZ = 0;

    // Add padding to ensure the entire mesh is inside the voxel grid
    // This prevents clipping of the mesh at the boundaries
    double padding = voxelSize;
    
    if (gridAligned) {
        std::cout << "\nUsing grid-aligned coordinates" << std::endl;
        
        // For grid-aligned mode, calculate offset in voxel units from global origin (0,0,0)
        offsetX = static_cast<int>(std::floor(minB.x / voxelSize));
        offsetY = static_cast<int>(std::floor(minB.y / voxelSize));
        offsetZ = static_cast<int>(std::floor(minB.z / voxelSize));

        // For grid-aligned mode, snap to voxel grid multiples
        minB.x = offsetX * voxelSize;
        minB.y = offsetY * voxelSize;
        minB.z = offsetZ * voxelSize;
        
        // Add padding to max bounds
        maxB.x = std::ceil((maxB.x + padding) / voxelSize) * voxelSize;
        maxB.y = std::ceil((maxB.y + padding) / voxelSize) * voxelSize;
        maxB.z = std::ceil((maxB.z + padding) / voxelSize) * voxelSize;
        
        std::cout << "  Adjusted bounds to align with voxel grid:" << std::endl;
        std::cout << "  X: " << minB.x << " to " << maxB.x << std::endl;
        std::cout << "  Y: " << minB.y << " to " << maxB.y << std::endl;
        std::cout << "  Z: " << minB.z << " to " << maxB.z << std::endl;

        std::cout << "Voxel grid offset from origin: " << offsetX << ", "
            << offsetY << ", " << offsetZ << " (in voxel units)" << std::endl;
    } else {
        // Add padding to ensure the entire mesh is inside the grid
        minB = minB - Vector3(padding, padding, padding);
        maxB = maxB + Vector3(padding, padding, padding);
    }

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
        std::queue<VoxelIdx> queue;
        
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
            VoxelIdx current = queue.front();
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
    std::string formatSuffix = denseOutput ? "-dense" : "-sparse";
    std::string outputFilename = baseFilename + formatSuffix + ".txt";
    
    // Call the metadata saving function before writing the actual voxel data
    saveVoxelMetadata(baseFilename, voxelSize, Nx, Ny, Nz, minB, maxB, offsetX, offsetY, offsetZ);

    auto writeStartTime = std::chrono::high_resolution_clock::now();

    // Optimized file writing with buffering for better performance
    // By using buffered I/O instead of writing each voxel individually, we minimize
    // the number of system calls and disk operations, which are typically much slower
    // than in-memory operations. This is especially important when dealing with
    // millions of voxels, as each system call has overhead.
    
    if (denseOutput) {
        // Write dense voxel format (full 3D array)
        std::cout << "Writing dense voxel format to " << outputFilename << std::endl;
        std::ofstream out(outputFilename);
        
        // Set larger buffer size (8MB) for better performance
        const size_t BUFFER_SIZE = 8 * 1024 * 1024;
        char* buffer = new char[BUFFER_SIZE];
        out.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);
        out.sync_with_stdio(false);
        
        // Write each voxel as 1 (occupied) or 0 (empty)
        int totalVoxels = Nx * Ny * Nz;
        int processedVoxels = 0;
        int lastPercentage = -1;
        
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                // Build each row as a string for better I/O performance
                std::string row;
                for (int i = 0; i < Nx; ++i) {
                    // Write 1 for occupied, 0 for empty
                    row += voxels[index(i, j, k)] ? "1" : "0";
                    if (i < Nx - 1) row += " "; // Space between values except last
                }
                out << row << "\n";
                
                // Update progress every 1%
                processedVoxels += Nx;
                int currentPercentage = (processedVoxels * 100) / totalVoxels;
                if (currentPercentage > lastPercentage) {
                    std::cout << "\rWriting voxels: " << currentPercentage << "% complete" << std::flush;
                    lastPercentage = currentPercentage;
                }
            }
        }
        
        // Clean up resources
        out.close();
        delete[] buffer;
    } else {
        // Write sparse voxel format (original implementation with i,j,k coordinates)
        std::cout << "Writing sparse voxel format to " << outputFilename << std::endl;
        std::ofstream out(outputFilename);
        
        // Set larger buffer size (8MB) to reduce system calls
        // This improves performance by batching write operations
        const size_t BUFFER_SIZE = 8 * 1024 * 1024;
        char* buffer = new char[BUFFER_SIZE];
        out.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);
        
        // Disable synchronization with C standard streams for better performance
        // This prevents unnecessary flushing between C and C++ I/O operations,
        // reducing overhead and improving file writing speed significantly
        out.sync_with_stdio(false);
        
        // Use a string stream for efficient string building before writing to file
        // A stringstream is an in-memory stream that allows for efficient string manipulation
        // without the overhead of file I/O operations. By accumulating output in memory first,
        // we can reduce the number of actual disk write operations, which are significantly
        // slower. Unlike writing directly to ofstream for each voxel, this approach batches
        // multiple voxels together before flushing to disk, greatly improving performance
        // when dealing with millions of voxels.
        std::stringstream ss;
        ss.sync_with_stdio(false);
        
        // Process in chunks to avoid excessive memory usage
        // This allows handling very large models without running out of memory
        const int CHUNK_SIZE = 1000000; // Process 1M voxels at a time
        std::vector<VoxelIdx> voxels_ijk;
        voxels_ijk.reserve(CHUNK_SIZE);
        
        // First pass: collect voxel centers in batches and write to file
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
                        voxels_ijk.push_back({i, j, k});
                        
                        // Write chunk if buffer is full to manage memory usage
                        if (voxels_ijk.size() >= CHUNK_SIZE) {
                            for (const auto& c : voxels_ijk) {
                                ss << c.i << "," << c.j << "," << c.k << "\n";
                            }
                            out << ss.str();  // Write the accumulated voxel data from stringstream to the output file
                            ss.str("");  // Clear the stringstream buffer to free memory and prepare for new data
                            voxels_ijk.clear();  // Clear the centers vector to free memory and prepare for next batch of voxels
                        }
                    }
                }
            }
        }
        
        // Write any remaining centers at the end
        if (!voxels_ijk.empty()) {
            for (const auto& c : voxels_ijk) {
                ss << c.i << "," << c.j << "," << c.k << "\n";
            }
            out << ss.str();
        }
        
        // Clean up resources
        out.close();
        delete[] buffer;
    }

    // Report file writing performance
    auto writeEndTime = std::chrono::high_resolution_clock::now();
    auto writeDuration = std::chrono::duration_cast<std::chrono::milliseconds>(writeEndTime - writeStartTime);
    std::cout << "\nVoxel centers saved to " << outputFilename << " in " 
              << writeDuration.count() / 1000.0 << " seconds\n";

    return 0;
}
