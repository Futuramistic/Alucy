#pragma once

// maximum number of vertices and triangles
#define MAXV 1000000
#define MAXT 1000000
#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <streambuf>

typedef int OrTri;
typedef int tIdx;

inline OrTri makeOrTri(tIdx t, int version) { return (t << 3) | version; };
inline tIdx idx(OrTri ot) { return ot >> 3; };
inline int ver(OrTri ot) { return ot & 0b111; };
inline OrTri enext(OrTri ot) {
	int v = ver(ot);  return makeOrTri(idx(ot),
	                           v < 3 ? (v + 1) % 3 : 3 + ((v - 1) % 3)) ; };
inline OrTri sym(OrTri ot) { int v = ver(ot); return v < 3 ? ot + 3 : ot - 3; };


class myObjType {
	int vcount = 0;

	double vlist[MAXV][3]; // vertices list
	double vlooplist[MAXV][3];
	int tlist[MAXT][3];      // triangle list
	int tlooplist[MAXT][3];
	int fnlist[MAXT][3];     // fnext list for future (not this assignment)
	double nlist[MAXT][3];   // storing triangle normals
	double vnlist[MAXV][3];  // storing vertex normals 
	

	
	std::vector<std::pair<int, int>> boundrylist; //storing boundry edges
	std::vector<std::set<int>> clist; //storing component lists

	std::set<int> neighbours[MAXV]; //storing vertex neighbours
	std::set<int> facelist[MAXV]; //storing facelist for each vertex
	std::set<int> boundryVertices;
	int collapseList[MAXV]; //vertex to collapse to
	double edgeCost[MAXV]; //cost of collapsing vertex

	double lmax[3];          // the maximum coordinates of x,y,z
	double lmin[3];          // the minimum coordinates of x,y,z
	int statMinAngle[18]; // each bucket is  degrees has a 10 degree range from 0 to 180 degree
	int statMaxAngle[18];
	double maxAngle = 0;
	double minAngle = 180;

	bool orientable = true;
	bool computedTriangleOrientation = false;
	bool boundry = false;

public:
	myObjType() { vcount = 0; tcount = 0; };
	int tcount=0;
	void readFile(char* filename);  // assumming file contains a manifold
	void loadSTL(char* filename);
	void orientTriangles();
	void writeFile(char* filename);
	void PrintInfo();
	void draw();
	void drawGouraud();
	void computeStat();
	void computeInfo();
	void computeTrianglesNormals();
	void computeVertexNormals();
	void computeBoundryEdges();
	std::vector<double> computeTriangleNormal(double vertex1[3], double vertex2[3], double vertex3[3]);
	void mergeTriangles(int i, int odd[3]);
	int computeOdd(int triangle,int j, std::map<std::pair<int, int>, int> &odds);
	void computeEven(int triangle);
	void loopSubdivide();
	void deleteVertex(int vertex);
	void deleteTriangle(int triangle);
	bool hasVertex(int vertex, int triangle);
	void computeAngles();
	void computeTriangleNormal(int i);
	float parseFloat(std::ifstream& s);
	void parseNormal(int triangle, std::ifstream& s);
	void parsePoint(int triangle, int vertexPos, std::ifstream& s);
	
	int org(OrTri ot);
	int dest(OrTri ot);
	void findFNext();
	void displayBoundries();
	void computeComponents();
	void getNeighbours();
	bool toggleBoundry();
	void simplifyMesh(int faceCount);
	void computeEdgeCost(int vertex);
	double computeEdgeCollapseCost(int vertex, int neighbour);
	void collapse(int vertex, int neighbour);
	void printNeighbours();
	void replace(int triangle,int vertex, int neighbour);
};


