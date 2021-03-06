#include "mesh.h"

#ifdef _WIN32
#include <Windows.h>
#include "GL\glut.h"
#define M_PI 3.141592654
#elif __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/GLUT.h>
#endif

#include "math.h"
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "mesh.h"
#include <map>
#include <queue>
#include <iomanip>
#include <set>
#include <queue>   

using namespace std;

/**
 ** Draws mesh onto the viewport 
 **/
void myObjType::draw() {
	glEnable(GL_LIGHTING);

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	glPushMatrix();
	double longestSide = 0.0;
	for (int i = 0; i < 3; i++)
		if ((lmax[i] - lmin[i]) > longestSide)
			longestSide = (lmax[i] - lmin[i]);
	glScalef(4.0 / longestSide, 4.0 / longestSide, 4.0 / longestSide);
	glTranslated(-(lmin[0] + lmax[0]) / 2.0, -(lmin[1] + lmax[1]) / 2.0, -(lmin[2] + lmax[2]) / 2.0);
	for (int i = 1; i <= tcount; i++)
	{
		glBegin(GL_POLYGON);
		glNormal3dv(nlist[i]);
		for (int j = 0; j < 3; j++) {
			if (Gouraud) {
				glNormal3dv(vnlist[tlist[i][j]]);
			}
			glVertex3dv(vlist[tlist[i][j]]);
		}
		glEnd();
	}	
	displayBoundries();
	glDisable(GL_LIGHTING);
	glPopMatrix();
}

/**
 ** Draws boundary edges on the mesh
 **/
void myObjType::displayBoundries() {
	if (boundry){
		float mat_ambient[] = { 1.0f, 0.0f, 0.0f, 1.0f };
		float mat_diffuse[] = { 1.0f, 0.0f, 0.0f, 1.0f };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glBegin(GL_LINES);
		for (int i = 0; i < boundrylist.size(); ++i) {
				glVertex3dv(vlist[boundrylist.at(i).first]);
				glVertex3dv(vlist[boundrylist.at(i).second]);
		}
		glEnd();
	}
}

/**
 ** Writes mesh into an .obj file
 ** param:
 ** - filename [*char] - name of a file to write to finished with '.obj'
 **/
void myObjType::writeFile(char* filename)
{
	std::ofstream outfile;
	outfile.open(filename, std::ios::out);
	std::string strbuff;
	for (int i = 1; i <= vcount; ++i) {
		outfile << "v " << vlist[i][0] << " " << vlist[i][1] << " " << vlist[i][2] << std::endl;
	}
	outfile << std::endl;
	for (int i = 1; i <= tcount; ++i) {
		outfile << "f " << tlist[i][0] << " " << tlist[i][1] << " " << tlist[i][2] << std::endl;
	}
	outfile.close();
}

/**
 ** Read mesh from an .obj file
 ** param:
 ** - filename [*char] - name of a file to read from finished with '.obj'
 **/
void myObjType::readFile(char* filename)
{
	cout << "Opening " << filename << endl;
	ifstream inFile;
	inFile.open(filename);
	if (!inFile.is_open()) {
		cout << "We cannot find your file " << filename << endl;
		exit(1);
	}

	string line;
	int i, j;
	bool firstVertex = 1;
	double currCood;

	while (getline(inFile, line))
	{
		if ((line[0] == 'v' || line[0] == 'f') && line[1] == ' ')
		{
			if (line[0] == 'v')
			{
				vcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) { // k is 0,1,2 for x,y,z
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ') j++;
					currCood = vlist[vcount][k] = atof(line.substr(i, j - i).c_str());
					if (firstVertex) 
						lmin[k] = lmax[k] = currCood;
					else {
						if (lmin[k] > currCood)
							lmin[k] = currCood;
						if (lmax[k] < currCood)
							lmax[k] = currCood;
					}
					i = j;
				}

				firstVertex = 0;
			}
			if (line[0] == 'f')
			{
				tcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) {
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ' && linec[j] != '\\') j++;
					tlist[tcount][k] = atof(line.substr(i, j - i).c_str());
					i = j;
					fnlist[tcount][k] = 0;
					while (linec[j] != ' ') j++;

				}

			}


		}
	}
	computeInfo();
}

/**
 ** Compute mesh information.
 ** Information include:
 ** normals, boundaries, components, angles, fnext table
 **/
void myObjType::computeInfo(){
	getNeighbours();
	computeTrianglesNormals();
	findFNext();
	computeComponents();
	computeBoundryEdges();
	computeVertexNormals();
	cout << "No. of vertices: " << vcount << endl;
	cout << "No. of triangles: " << tcount << endl;
	computeStat();
}

/**
 ** Find origin of oriTri representation
 ** param:
 ** - ot [OrTri]: oriTri number
 **/
int myObjType::org(OrTri ot) {
	int v = ver(ot);
	int id = idx(ot);
	if (v < 3) {
		return tlist[id][v];
	}
	else {
		return tlist[id][(v+1)%3];
	}
}

/**
 ** Find destination of oriTri representation
 ** param:
 ** - ot [OrTri]: oriTri number
 **/
int myObjType::dest(OrTri ot) {
	int v = ver(ot);
	int id = idx(ot);
	if (v < 3) {
		return tlist[id][(v+1)%3];
	}
	else {
		return tlist[id][v % 3];
	}

}

/**
 ** Compute fnext table
 **/
void myObjType::findFNext(){
	for (int i = 1; i <=tcount; ++i) {
		for (int j = 0; j < 3; ++j) {
			fnlist[i][j] = 0;
		}
	}
	std::map<std::pair<int, int>, int> hashMap;
	for (int i = 1; i <= tcount; ++i){
		for (int j = 0; j < 3; ++j) {
			int oriTri = makeOrTri(i, j);
			int start = org(oriTri), end = dest(oriTri);
			std::pair <int, int>key1(start, end);	
			std::map<pair <int, int>, int>::iterator it1 = hashMap.find(key1);
			if (it1 != hashMap.end()) {
				int found_oriTri = hashMap.at(key1);
				int id = idx(found_oriTri);
				int v = ver(found_oriTri);
				fnlist[i][j]= found_oriTri;
				fnlist[id][v%3] = oriTri;
				hashMap.erase(it1);
			}
			else{
				std::pair <int, int>key2(end, start);
				std::map<pair <int, int>, int>::iterator it2 = hashMap.find(key2);
				if (it2 != hashMap.end()) {
					int found_oriTri = hashMap.at(key2);
					int id = idx(found_oriTri);
					int v = ver(found_oriTri);
					fnlist[i][j] = found_oriTri;
					fnlist[id][v%3] = sym(oriTri);
					hashMap.erase(it2);
				}
				else {
					hashMap.emplace(key1, oriTri);
				}
			}
		}
	}

	for (int i = 1; i <= tcount; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (fnlist[i][j] == 0) {
				fnlist[i][j] = makeOrTri(i, j);
			}
		}
	}
	
}

/**
 ** Compute all triangles normals
 **/
void myObjType::computeTrianglesNormals() {
	for (int i = 1; i <= tcount; ++i) {
		computeTriangleNormal(i);
	}
}

/**
 ** Compute triangle normal
 ** param:
 ** - i [int]: triangle index
 **/
void myObjType::computeTriangleNormal(int i) {
	double ax = vlist[tlist[i][1]][0] - vlist[tlist[i][0]][0];
	double ay = vlist[tlist[i][1]][1] - vlist[tlist[i][0]][1];
	double az = vlist[tlist[i][1]][2] - vlist[tlist[i][0]][2];

	double bx = vlist[tlist[i][2]][0] - vlist[tlist[i][0]][0];
	double by = vlist[tlist[i][2]][1] - vlist[tlist[i][0]][1];
	double bz = vlist[tlist[i][2]][2] - vlist[tlist[i][0]][2];

	double nx = (ay * bz) - (az * by);
	double ny = (az * bx) - (ax * bz);
	double nz = (ax * by) - (ay * bx);

	double length = sqrt(nx * nx + ny * ny + nz * nz);
	if (length > 0.0) {
		nx /= length;
		ny /= length;
		nz /= length;
	}
	nlist[i][0] = nx;
	nlist[i][1] = ny;
	nlist[i][2] = nz;
}

/**
 ** Compute triangle normal
 ** param:
 ** - vertex1, vertex2, vertex3 [double[3]] - vertices
 **/
std::vector<double> myObjType::computeTriangleNormal(double vertex1[3], double vertex2[3], double vertex3[3]) {
	double ax = vertex2[0] - vertex1[0];
	double ay = vertex2[1] - vertex1[1];
	double az = vertex2[2] - vertex1[2];

	double bx = vertex3[0] - vertex1[0];
	double by = vertex3[1] - vertex1[1];
	double bz = vertex3[2] - vertex1[2];

	double nx = (ay * bz) - (az * by);
	double ny = (az * bx) - (ax * bz);
	double nz = (ax * by) - (ay * bx);

	double length = sqrt(nx * nx + ny * ny + nz * nz);
	if (length > 0.0) {
		nx /= length;
		ny /= length;
		nz /= length;
	}
	return std::vector<double>{nx, ny, nz};
}


/**
 ** Compute mesh angles statistics
 **/
void myObjType::computeAngles() {
	for (int i = 0; i < 18; ++i) {
		statMinAngle[i] = 0;
		statMaxAngle[i] = 0;
	}
	maxAngle = 0;
	minAngle = 180;
	for (int i = 1; i <= tcount; ++i) {
		int minBucket=18;
		int maxBucket = 0;
		for (int j = 0; j < 3; ++j) {
			int first;
			int second;
			switch (j) {
				case 0: {
					first = 1;
					second = 2;
					break;
				}
				case 1: {
					first = 0;
					second = 2;
					break;
				}
				case 2: {
					first = 1;
					second = 0;
					break;
				}

			}
			double ax = vlist[tlist[i][first]][0] - vlist[tlist[i][j]][0];
			double ay = vlist[tlist[i][first]][1] - vlist[tlist[i][j]][1];
			double az = vlist[tlist[i][first]][2] - vlist[tlist[i][j]][2];
			double alength = sqrt(ax * ax + ay * ay + az * az);

			double bx = vlist[tlist[i][second]][0] - vlist[tlist[i][j]][0];
			double by = vlist[tlist[i][second]][1] - vlist[tlist[i][j]][1];
			double bz = vlist[tlist[i][second]][2] - vlist[tlist[i][j]][2];
			double blength = sqrt(bx * bx + by * by + bz * bz);

			float alpha = acos((ax * bx + ay * by + az * bz) / (alength * blength)) * 180 / M_PI;
			int bucket = floor(alpha / 10);
			if (bucket < minBucket) {
				minBucket = bucket;
			}
			if (bucket > maxBucket) {
				maxBucket = bucket;
			}
			if (alpha > maxAngle) {
				maxAngle = alpha;
			}
			if (alpha < minAngle) {
				minAngle = alpha;
			}
		}
		++statMinAngle[minBucket];
		++statMaxAngle[maxBucket];
	}
}

/**
 ** Compute mesh components
 **/
void myObjType::computeComponents() {
	clist.clear();
	std::set<int> triangles;
	for (int i = 1; i <= tcount; ++i) {
		triangles.insert(i);
	}
	while (!triangles.empty()) {
		std::priority_queue<int> queue;
		std::set<int> compontentTriangles;
		int first = *triangles.begin();
		queue.push(first);
		compontentTriangles.insert(first);
		while (!queue.empty()) {
			int index = queue.top();
			queue.pop();
			for (int j = 0; j < 3; ++j) {
				if (compontentTriangles.find(idx(fnlist[index][j])) == compontentTriangles.end() && 
					triangles.find(idx(fnlist[index][j])) != triangles.end()){
					queue.push(idx(fnlist[index][j]));
					compontentTriangles.insert(idx(fnlist[index][j]));
				}
			}
			triangles.erase(index);
		}
		clist.push_back(compontentTriangles);
	}
}

/**
 ** Compute and display mesh statistics
 **/
void myObjType::computeStat()
{
	int i;
	computeAngles();
    cout << "Min. angle = " << minAngle << endl;
    cout << "Max. angle = " << maxAngle << endl;
	
	cout << "Statistics for Maximum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMaxAngle[i] << " ";
	cout << endl;
	cout << "Statistics for Minimum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMinAngle[i] << " ";
	cout << endl;
	cout << "Components: " << clist.size()<<endl;
	cout << "Boundry edges: " << boundrylist.size() << endl;
}

/**
 ** Orient triangles
 **/
void myObjType::orientTriangles(){
	if (!orientable) {
		cout << "ERROR: MESH CANNOT BE ORIENTED" << endl;
		return;
	}
	set<int> flipped;
	for (int i = 1; i <= tcount; ++i) {
		//Get this triangle info
		int ori = makeOrTri(i, 0);
		int start = org(ori);
		int end = dest(ori);
		//Get neighbour triangle info
		int next = fnlist[i][0];
		int nextIdx = idx(next);
		int nextStart = org(next);
		int nextEnd = dest(next);
		if (nextEnd == end || nextStart == start && nextIdx != i)
		{	
			if (flipped.find(nextIdx) != flipped.end()){
				orientable = false;
				cout <<"ERROR: MESH CANNOT BE ORIENTED"<<endl;
				return;
			}
			if (nextIdx < i){
				int temp = tlist[nextIdx][0];
				tlist[nextIdx][0] = tlist[nextIdx][1];
				tlist[nextIdx][1] = temp;
				for (int j = 0; j < 3; ++j) {
					fnlist[nextIdx][j] = sym(fnlist[nextIdx][j]);
				}
				fnlist[i][0] = sym(fnlist[i][0]);
				computeTriangleNormal(nextIdx);
				flipped.insert(nextIdx);
			}
		}
	}
	findFNext();
	computeVertexNormals();
}

/**
 ** Make the boundary visible or not
 **/
bool myObjType::toggleBoundry() {
	boundry = !boundry;
	return boundry;
}

/**
 ** Compute normals for each vertex
 **/
void myObjType::computeVertexNormals() {
	for (int i = 1; i <= vcount; ++i) {
		for (int j = 0; j < 3; ++j) {
			vnlist[i][j] = 0;
		}
	}
	for (int i = 1; i <= tcount; ++i) {
		for (int j = 0; j < 3; ++j){
			int triangle = i;
			int vertex = org(makeOrTri(i, j));
			int next = enext(sym(fnlist[triangle][j]));
			set<int> triangles;
			while (triangles.find(triangle) == triangles.end()) {
					vnlist[vertex][0] += nlist[triangle][0];
					vnlist[vertex][1] += nlist[triangle][1];
					vnlist[vertex][2] += nlist[triangle][2];
					triangles.insert(triangle);
					triangle = idx(next);
					next = enext(sym(fnlist[triangle][j]));
			}
			if (triangles.size() > 0) {
					vnlist[vertex][0]/=triangles.size();
					vnlist[vertex][1]/=triangles.size();
					vnlist[vertex][2]/=triangles.size();
			}
			double length = sqrt(vnlist[vertex][0] * vnlist[vertex][0] + vnlist[vertex][1] * vnlist[vertex][1] + vnlist[vertex][2] * vnlist[vertex][2]);
			if (length > 0.0) {
					vnlist[vertex][0] /= length;
					vnlist[vertex][1] /= length;
					vnlist[vertex][2] /= length;
			}
		}
	}
}

/**
 ** Find boudary edges
 **/
void myObjType::computeBoundryEdges(){
	boundrylist.clear();
	boundryVertices.clear();
	for (int i = 1; i <= tcount; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (idx(fnlist[i][j]) == i) {
				int orTri = makeOrTri(i, j);
				boundrylist.push_back(std::pair<int,int>(org(orTri),dest(orTri)));
				boundryVertices.insert(org(orTri));
				boundryVertices.insert(dest(orTri));
			}
		}
	}
}

/**
 ** Compute neighbours for each vertex
 **/
void myObjType::getNeighbours() {
	for (int i = 1; i <= vcount; ++i) {
		neighbours[i].clear();
		facelist[i].clear();
	}
	for (int i = 1; i <= tcount; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				if (tlist[i][j]!=tlist[i][k]){
					neighbours[tlist[i][j]].insert(tlist[i][k]);
				}
			}
			facelist[tlist[i][j]].insert(i);
		}
	}
}

/**
 ** Simplify mesh to desired number of faces
 ** param:
 ** -faceCount [int] - number of faces to simplify to
 **/
void myObjType::simplifyMesh(int faceCount) {
	while (tcount > faceCount) {
		getNeighbours();
		computeTrianglesNormals();
		for (int i = 1; i <= vcount; ++i) {
			computeEdgeCost(i);
		}
		int vertex = 1;
		for (int i = 1; i <= vcount; ++i) {
			if (edgeCost[i] < edgeCost[vertex]) {
				vertex = i;
			}
		}
		collapse(vertex, collapseList[vertex]);
	}
	orientable = 1;
	computeInfo();
}

/**
 ** Delete vertex from mesh
 ** param:
 ** - vertex [int]: vertex index
 **/
void myObjType::deleteVertex(int vertex){
	if (vertex > 0) {
		//Move all indexes by one if bigger than vertex
		for (int i = 1;i<=tcount;++i){
			for(int j = 0;j<3;++j){
				if (tlist[i][j]>vertex){
					--tlist[i][j];
				}
			}
		}
		//Move all vertices
		for(int i=vertex;i<vcount;i++){
			for (int j = 0; j < 3; ++j){
				vlist[i][j] = vlist[i + 1][j];
			}
		}
		--vcount;
	}
	return;
}

/**
 ** Check if triangle has vertex
 ** param:
 ** -vertex [int]: vertex index to check for
 ** -triangle [int]: triangle index to check in
 **/
bool myObjType::hasVertex(int vertex, int triangle) {
	return (tlist[triangle][0] == vertex || tlist[triangle][1] == vertex || tlist[triangle][2] == vertex);
}

/**
 ** Delete triangle from mesh
 ** -triangle [int]: triangle index
 **/
void myObjType::deleteTriangle(int triangle) {
	//Move all the triangles
	if (triangle > 0) {
		for (int i = triangle; i <= tcount - 1; ++i) {
			for (int j = 0; j < 3; ++j) {
				tlist[i][j] = tlist[i + 1][j];
			}
		}
		--tcount;
	}
}

/**
 ** Half-collapse edge
 ** param:
 ** -vertex [int]: vertex to be collapsed
 ** -neighbour [int]: vertex to collapse into
 **/
void myObjType::collapse(int vertex, int neighbour) {
	if (neighbour == -1) {
		deleteVertex(vertex);
		return;
	}
	std::vector<int> deletedFaces;
	for (set<int>::iterator it = facelist[vertex].begin(); it!= facelist[vertex].end(); ++it) {
		int triangle = *it;
		if (hasVertex(neighbour, triangle)) {deletedFaces.push_back(triangle);}
		else{replace(triangle,vertex,neighbour);}
	}
	for (int i = 0; i<deletedFaces.size(); ++i) {
		int triangle = deletedFaces.at(i);
		deleteTriangle(triangle);
		for (int j = i+1; j < deletedFaces.size(); ++j) {
			if (triangle < deletedFaces.at(j)){
				--deletedFaces[j];
			}
		}
	}
	deleteVertex(vertex);
}

/**
 ** Replace vertex with neighbour in triangle
 ** param:
 ** -triangle [int]: triangle index
 ** -vertex [int]: vertex to be replaced
 ** -neighbour [int]: vertex to replace with
 **/
void myObjType::replace(int triangle,int vertex, int neighbour) {
	if(tlist[triangle][0]==vertex){
		tlist[triangle][0] = neighbour;
	}
	else if (tlist[triangle][1]==vertex){
		tlist[triangle][1] = neighbour;
	}
	else{
		tlist[triangle][2] = neighbour;
	}
}

/**
 ** Compute smallest edge cost for vertex
 ** param:
 ** -vertex [int]: vertex to compute the cost for
 **/
void myObjType::computeEdgeCost(int vertex) {
	if (neighbours[vertex].size() == 0) {
		collapseList[vertex] = -1;
		edgeCost[vertex] = -0.1;
		return;
	}
	collapseList[vertex] = -1;
	edgeCost[vertex] = 100000000;
	for (set<int>::iterator i = neighbours[vertex].begin(); i != neighbours[vertex].end(); ++i) {
		double cost = computeEdgeCollapseCost(vertex, *i);
		if (cost < edgeCost[vertex]) {
			collapseList[vertex] = *i;
			edgeCost[vertex] = cost;
		}
	}
}

/**
 ** Compute edge collapse cost
 ** param:
 ** -vertex, neighbour [int]: edge to compute for
 **/
double myObjType::computeEdgeCollapseCost(int vertex, int neighbour) {
	double edgeLength = sqrt((vlist[vertex][0]-vlist[neighbour][0])*(vlist[vertex][0]-vlist[neighbour][0]) + 
							 (vlist[vertex][1]-vlist[neighbour][1])*(vlist[vertex][1]-vlist[neighbour][1]) + 
							 (vlist[vertex][2]-vlist[neighbour][2])*(vlist[vertex][2]-vlist[neighbour][2]));
	double curvature = 0.0;

	std::vector<int> sides;
	for (set<int>::iterator it = facelist[vertex].begin(); it != facelist[vertex].end(); ++it) {
		int triangle = *it;
		if (hasVertex(neighbour,triangle)){ sides.push_back(triangle); }
	}
	for (set<int>::iterator it = facelist[vertex].begin(); it != facelist[vertex].end(); ++it) {
		double minCurv = 1.0;
		for (int j = 0; j < sides.size(); j++) {
			int triangle = sides.at(j);
			double dotProd = (nlist[*it][0]*nlist[triangle][0])+(nlist[*it][1]*nlist[triangle][1]) + (nlist[*it][2] * nlist[triangle][2]);
			if (minCurv > (1 - dotProd) / 2.0){ minCurv = (1 - dotProd) / 2.0;}
		}
		if (curvature<minCurv){curvature = minCurv;}
	}
	return (edgeLength*curvature)/(std::pow(facelist[vertex].size(),2));
}

/**
 ** Copy coordinates of even vertices
 ** param:
 ** -traingle [int]: index of triangle to copy the vertices from
 **/
void myObjType::copyEven(int triangle) {
	//Don't move the previous vertices
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			vlooplist[tlist[triangle][i]][j] = vlist[tlist[triangle][i]][j];
		}
	}
}

/**
 ** Compute even vertices based on loop subdivison algorithm
 ** param:
 ** -traingle [int]: index of triangle to compute the vertices for
 **/
void myObjType::computeEven(int triangle) {
	for (int j = 0; j < 3; ++j) {
		if (vlooplist[tlist[triangle][j]][0] == 0) {
			bool crease = false;
			//check if boundary vertex
			int vertex = tlist[triangle][j];
			double angle = 1.0;
			int vertices[2];
			std::vector<int> sides;
			std::vector<int> nei;
			for (set<int>::iterator it = neighbours[vertex].begin(); it != neighbours[vertex].end(); ++it) {
				int neighbour = *it;
				for (set<int>::iterator it = facelist[neighbour].begin(); it != facelist[neighbour].end(); ++it) {
					if (hasVertex(vertex, *it) && *it!=triangle) {
						sides.push_back(triangle);
						nei.push_back(neighbour);
					}
				}
			}
			for (set<int>::iterator it = facelist[vertex].begin(); it != facelist[vertex].end(); ++it) {
				double minCos = 1.0;
				for (int j = 0; j < sides.size(); j++) {
					int side = sides.at(j);
					double dotProd = (nlist[*it][0] * nlist[side][0]) + (nlist[*it][1] * nlist[side][1]) + (nlist[*it][2] * nlist[side][2]);
					double cos =dotProd/(sqrt(nlist[*it][0] * nlist[*it][0] + nlist[*it][1] * nlist[*it][1] + nlist[*it][2] * nlist[*it][2]) *
							   sqrt(nlist[side][0] * nlist[side][0] + nlist[side][1] * nlist[side][1] + nlist[side][2] * nlist[side][2]));
					if (abs(minCos) > abs(cos)) {
						minCos = cos;
						vertices[0] = nei.at(j);
					};
				}
				if (abs(angle) > abs(minCos)) { angle = minCos; }
			}
			if (abs(angle) < 0.99) {
				crease = true;
			}
			if (!crease) {
				int k = neighbours[tlist[triangle][j]].size();
				double beta = 1 / k * (5 / 8 - std::pow(2, (3 / 8 + 1 / 4 * cos(2 * M_PI / k))));
				double x = 0;
				double y = 0;
				double z = 0;
				for (set<int>::iterator it = neighbours[tlist[triangle][j]].begin(); it != neighbours[tlist[triangle][j]].end(); ++it) {
					x += vlist[*it][0];
					y += vlist[*it][1];
					z += vlist[*it][2];
				}
				double vx = vlist[tlist[triangle][j]][0];
				double vy = vlist[tlist[triangle][j]][1];
				double vz = vlist[tlist[triangle][j]][2];
				vlooplist[tlist[triangle][j]][0] = vx * (1.0 - (k * beta)) + x * beta;
				vlooplist[tlist[triangle][j]][1] = vy * (1.0 - (k * beta)) + y * beta;
				vlooplist[tlist[triangle][j]][2] = vz * (1.0 - (k * beta)) + z * beta;
			}
			else {
				int a = vertices[0];
				int b;
				for (int i = 0; i < 3; ++i) {
					if (tlist[triangle][i] != a && tlist[triangle][i] != vertex) {
						b = tlist[triangle][i];
					}
				}
				double vx = vlist[tlist[triangle][j]][0];
				double vy = vlist[tlist[triangle][j]][1];
				double vz = vlist[tlist[triangle][j]][2];
				vlooplist[tlist[triangle][j]][0] = 0.75 * vx + 0.125 * (vlist[a][0]+ vlist[b][0]);
				vlooplist[tlist[triangle][j]][1] = 0.75 * vy + 0.125 * (vlist[a][1] + vlist[b][1]);
				vlooplist[tlist[triangle][j]][2] = 0.75 * vz + 0.125 * (vlist[a][2] + vlist[b][2]);
			}
		}
	}
}

/**
 ** Compute the odd vertices based on loop algorithm
 ** param:
 ** -traingle [int]: index of triangle to copy the vertices from
 ** -odds [map<pair<int,int>,int>]: previously computed odd vertices
 **/
int myObjType::computeOdd(int triangle, int j, std::map<std::pair<int, int>, int> &odds) {
	if (idx(fnlist[triangle][j]) == triangle) {
			++vcount;
			int a = org(fnlist[triangle][j]);
			int b = dest(fnlist[triangle][j]);
			vlooplist[vcount][0] = 1.0 / 2.0 * (vlist[a][0] + vlist[b][0]);
			vlooplist[vcount][1] = 1.0 / 2.0 * (vlist[a][1] + vlist[b][1]);
			vlooplist[vcount][2] = 1.0 / 2.0 * (vlist[a][2] + vlist[b][2]);
			return vcount;
	}
	else{
			int a = org(fnlist[triangle][j]);
			int b = dest(fnlist[triangle][j]);
			int c, d;
			for (int i = 0; i < 3; ++i) {
				if (tlist[idx(fnlist[triangle][j])][i] != a && tlist[idx(fnlist[triangle][j])][i] != b) {
					c = tlist[idx(fnlist[triangle][j])][i];
				}
			}
			for (int i = 0; i < 3; ++i) {
				if (tlist[triangle][i] != a && tlist[triangle][i] != b) {
					d = tlist[triangle][i];
				}
			}
			std::pair<int, int> pair = { c,d };
			std::map<std::pair<int, int>, int>::iterator it = odds.find(pair);
			if (it!=odds.end()) {
				int value = odds.at(pair);
				return value;
			}
			else {
				pair = { d,c };
				std::map<std::pair<int, int>, int>::iterator it = odds.find(pair);
				if(it != odds.end()){
					int value = odds.at(pair);
					return value;
				}
				else {
					++vcount;
					vlooplist[vcount][0] = (3.0 / 8.0) * (vlist[a][0] + vlist[b][0]) + (1.0 / 8.0) * (vlist[c][0] + vlist[d][0]);
					vlooplist[vcount][1] = (3.0 / 8.0) * (vlist[a][1] + vlist[b][1]) + (1.0 / 8.0) * (vlist[c][1] + vlist[d][1]);
					vlooplist[vcount][2] = (3.0 / 8.0) * (vlist[a][2] + vlist[b][2]) + (1.0 / 8.0) * (vlist[c][2] + vlist[d][2]);
					odds.emplace(pair, vcount);
					return vcount;
				}	
			}
	}
}

/**
 ** Subdivide the mesh based on loop subdivision algorithm
 **/
void myObjType::loopSubdivide(){
	std::map<std::pair<int, int>, int> odds =  std::map<std::pair<int, int>, int>();
	getNeighbours();
	int odd[3];
	int triangles = tcount;
	for (int i = 1; i <= triangles; ++i) {
		for (int j = 0; j < 3; ++j) {
			odd[j]=computeOdd(i,j,odds);
		}
		copyEven(i);
		mergeTriangles(i,odd);
	}
	clearLoop();
	computeInfo();
}

/**
 ** Clear loop arrays
 **/
void myObjType::clearLoop() {
	for (int i = 1; i <= vcount; ++i) {
		vlist[i][0] = vlooplist[i][0];
		vlist[i][1] = vlooplist[i][1];
		vlist[i][2] = vlooplist[i][2];
		vlooplist[i][0] = 0;
		vlooplist[i][1] = 0;
		vlooplist[i][2] = 0;
	}
	for (int i = 1; i <= tcount; ++i) {
		tlist[i][0] = tlooplist[i][0];
		tlist[i][1] = tlooplist[i][1];
		tlist[i][2] = tlooplist[i][2];
		tlooplist[i][0] = 0;
		tlooplist[i][1] = 0;
		tlooplist[i][2] = 0;
	}
}

/**
 ** Make traingles for loop subdivision algorithm
 **/
void myObjType::mergeTriangles(int i, int odd[3]){
	int odd1 = odd[0];
	int odd2 = odd[1];
	int odd3 = odd[2];
	int even1 = tlist[i][0];
	int even2 = tlist[i][1];
	int even3 = tlist[i][2];
	tlooplist[i][0] = odd1;
	tlooplist[i][1] = odd2;
	tlooplist[i][2] = odd3;
	std::vector<double> normal = computeTriangleNormal(vlooplist[odd1], vlooplist[odd2], vlooplist[odd3]);
	double nx = normal.at(0);
	double ny = normal.at(1);
	double nz = normal.at(2);
	if (nx * nlist[i][0] + ny * nlist[i][1] + nz * nlist[i][2] < 0) {
		tlooplist[i][0] = odd2;
		tlooplist[i][1] = odd1;
		tlooplist[i][2] = odd3;
	}
	++tcount;
	tlooplist[tcount][0] = odd2;
	tlooplist[tcount][1] = odd1;
	tlooplist[tcount][2] = even2;
	normal = computeTriangleNormal(vlooplist[odd2], vlooplist[odd1], vlooplist[even2]);
	nx = normal.at(0);
	ny = normal.at(1);
	nz = normal.at(2);
	if (nx * nlist[i][0] + ny * nlist[i][1] + nz * nlist[i][2] < 0) {
		tlooplist[tcount][0] = odd1;
		tlooplist[tcount][1] = odd2;
		tlooplist[tcount][2] = even2;
	}
	++tcount;
	tlooplist[tcount][0] = odd1;
	tlooplist[tcount][1] = odd3;
	tlooplist[tcount][2] = even1;
	normal = computeTriangleNormal(vlooplist[odd1], vlooplist[odd3], vlooplist[even1]);
	nx = normal.at(0);
	ny = normal.at(1);
	nz = normal.at(2);
	if (nx * nlist[i][0] + ny * nlist[i][1] + nz * nlist[i][2] < 0) {
		tlooplist[tcount][0] = odd3;
		tlooplist[tcount][1] = odd1;
		tlooplist[tcount][2] = even1;
	}
	++tcount;
	tlooplist[tcount][0] = odd3;
	tlooplist[tcount][1] = odd2;
	tlooplist[tcount][2] = even3;
	normal = computeTriangleNormal(vlooplist[odd3], vlooplist[odd2], vlooplist[even3]);
	nx = normal.at(0);
	ny = normal.at(1);
	nz = normal.at(2);
	if (nx * nlist[i][0] + ny * nlist[i][1] + nz * nlist[i][2] < 0) {
		tlooplist[tcount][0] = odd2;
		tlooplist[tcount][1] = odd3;
		tlooplist[tcount][2] = even3;
	}
}

/**
 ** Load an ASCII STL file
 ** param:
 ** -fname [char*]: name of an stl file to load (should finish with .stl)
 **/
void myObjType::loadSTL(char* fname) {
	ifstream myFile(fname);
	if (!myFile) {
		cout << "We cannot find your file " << fname << endl;
		exit(1);
	}
	string buffer;
	vector<string> tokens;
	int maxNumTokens = 0;
	bool firstVertex = 1;
	int vertices[3];
	int triangleVertex = 0;
	vcount = 1;
	tcount = 1;
	while (!(myFile.eof() || myFile.fail()))
	{
		getline(myFile, buffer);
		istringstream line(buffer);
		int tokenCount = 0;
		while (!(line.eof() || line.fail())) {
			if (tokenCount >= maxNumTokens) {
				maxNumTokens = tokenCount + 1;
				tokens.resize(maxNumTokens);
			}
			line >> tokens[tokenCount];
			++tokenCount;
		}
		if (tokenCount > 0)
		{
			string& tok = tokens[0];
			if (tok.compare("vertex")== 0){
				double c[3];
				for (int i = 0; i < 3; ++i) {
					double currCood = static_cast<double> (atof(tokens[i + 1].c_str()));
					if (firstVertex)
						lmin[i] = lmax[i] = currCood;
					else {
						if (lmin[i] > currCood)
							lmin[i] = currCood;
						if (lmax[i] < currCood)
							lmax[i] = currCood;
					}
					c[i] = currCood;
				}
				bool found = false;
				int vertexFound;
				for (int i = 1; i < vcount; ++i) {
					if (vlist[i][0] == c[0] && vlist[i][1] == c[1] && vlist[i][2] == c[2]) {
						found = true;
						vertexFound = i;
					}
				}
				if (!found) {
					for (int i = 0; i < 3; ++i) {
						vlist[vcount][i] = c[i];
					}
					vertices[triangleVertex] = vcount;
					++vcount;
				}
				else {
					vertices[triangleVertex] = vertexFound;
				}
				firstVertex = 0;
				++triangleVertex;
			}
			else if (tok.compare("facet")==0){
			}
			else if (tok.compare("outer")==0 || tok.compare("endsolid") == 0) {
			}
			else if (tok.compare("endfacet") == 0){
				for (int i = 0; i < 3; ++i) {
					tlist[tcount][i] = vertices[i];
				}
				triangleVertex = 0;
				++tcount;
			}
		}
	}
	computeInfo();
	return;
}