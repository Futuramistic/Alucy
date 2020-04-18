#pragma once

// maximum number of vertices and triangles
#define MAXV 1000000
#define MAXT 1000000
#include <map>
#include <set>
#include <vector>

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
	int tcount = 0;

	double vlist[MAXV][3];   // vertices list
	int tlist[MAXT][3];      // triangle list
	int fnlist[MAXT][3];     // fnext list for future (not this assignment)
	double nlist[MAXT][3];   // storing triangle normals
	double vnlist[MAXV][3];
	std::vector<std::pair<int,int>> boundrylist;
	bool boundry = false;
	std::vector<std::set<int>> clist;
	double lmax[3];          // the maximum coordinates of x,y,z
	double lmin[3];          // the minimum coordinates of x,y,z

	int statMinAngle[18]; // each bucket is  degrees has a 10 degree range from 0 to 180 degree
	int statMaxAngle[18]; 
	float maxAngle = 0;
	float minAngle = 180;
	int components;
	bool orientable = true;
	bool computedTriangleOrientation = false;

	std::map<std::pair<int, int>, int> hashMap;

public:
	myObjType() { vcount = 0; tcount = 0; };
	void readFile(char* filename);  // assumming file contains a manifold
	void load3DS(char* filename);
	void orientTriangles();
	void writeFile(char* filename);  
	void OriTriPrint();
	void draw();  
	void drawGouraud();
    void computeStat();
	void computeTrianglesNormals();
	void computeVertexNormals();
	void computeBoundryEdges();
	void computeAngles();
	void computeTriangleNormal(int i);
	int org(OrTri ot);
	int dest(OrTri ot);
	void findFNext();
	void displayBoundries();
	void computeComponents();
	void toggleBoundry();
};


