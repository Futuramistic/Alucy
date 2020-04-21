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
// uncomment the following after you computed the normals
		glNormal3dv(nlist[i]);
		for (int j = 0; j < 3; j++) {
			glVertex3dv(vlist[tlist[i][j]]);
		}
		glEnd();
	}	
	displayBoundries();
	glDisable(GL_LIGHTING);
	glPopMatrix();
}

void myObjType::drawGouraud() {
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
		for (int j = 0; j < 3; j++) {
			glNormal3dv(vnlist[tlist[i][j]]);
			glVertex3dv(vlist[tlist[i][j]]);
		}
		glEnd();

	}
	displayBoundries();
	glDisable(GL_LIGHTING);
	glPopMatrix();
}

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

	// We suggest you to compute the normals here
	computeInfo();
}

void myObjType::computeInfo(){
	computeTrianglesNormals();
	findFNext();
	computeComponents();
	computeBoundryEdges();
	computeVertexNormals();
	cout << "No. of vertices: " << vcount << endl;
	cout << "No. of triangles: " << tcount << endl;
	computeStat();
}

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

void myObjType::findFNext(){
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
					fnlist[id][v%3] = oriTri;
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

void myObjType::computeTrianglesNormals() {
	for (int i = 1; i <= tcount; ++i) {
		computeTriangleNormal(i);
	}
}

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



void myObjType::computeAngles() {
	for (int i = 1; i <= tcount; ++i) {
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
			++statMaxAngle[bucket];
			++statMinAngle[bucket];
			if (alpha > maxAngle) {
				maxAngle = alpha;
			}
			if (alpha < minAngle) {
				minAngle = alpha;
			}
		}
	}
}

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

void myObjType::PrintInfo() {
	std::ofstream outfile;
	outfile.open("log.txt", std::ios::out);
	std::string strbuff;
	for (int i = 1; i <= tcount; ++i) {
		outfile << "Triangle: " << i<<": "<<endl;
		for (int j = 0; j < 3; ++j) {
			 outfile<<"Vertex: "<< tlist[i][j] << endl;
		}
		outfile << endl;
	}
	outfile.close();
}

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

void myObjType::orientTriangles(){
	if (!orientable) {
		cout << "Not orientable" << endl;
		return;
	}
	for (int i = 0; i < clist.size(); ++i) {
		std::set<int> component = clist[i];
		int start = *component.begin();
		std::priority_queue<int> queue;
		queue.push(start);
		while (!queue.empty() && orientable) {
			int index = queue.top();
			queue.pop();
			for (int j = 0; j < 3; ++j) {
				int next = idx(fnlist[index][j]);
				double dot = nlist[next][0] * nlist[index][0] + nlist[next][1] * nlist[index][1] + nlist[next][2] * nlist[index][2];
				if (dot < 0.0) {
					if (component.find(next) != component.end()) {
						queue.push(next);
						int temp = tlist[next][0];
						tlist[next][0] = tlist[next][1];
						tlist[next][1] = temp;
						computeTriangleNormal(next);
					}
					else {
						orientable = false;
						break;
					}
				}
			}
			component.erase(index);
		}
		if (!orientable) {
			cout << "Not orientable" << endl;
			return;
		}
	}
		computeVertexNormals();
	}

void myObjType::toggleBoundry() {
	boundry = !boundry;
}

void myObjType::computeVertexNormals() {
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
	computeInfo();
}

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

bool myObjType::hasVertex(int vertex, int triangle) {
	return (tlist[triangle][0] == vertex || tlist[triangle][1] == vertex || tlist[triangle][2] == vertex);
}

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
	return edgeLength*curvature;
}

void myObjType::computeEven(int triangle) {
	for (int j = 0; j < 3; ++j) {
		if (vlooplist[tlist[triangle][j]][0] == 0) {
			bool crease = false;
			//check if boundary vertex
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
				double vx = vlist[tlist[triangle][j]][0];
				double vy = vlist[tlist[triangle][j]][1];
				double vz = vlist[tlist[triangle][j]][2];
				vlooplist[tlist[triangle][j]][0] = 0.75 * vx;
				vlooplist[tlist[triangle][j]][1] = 0.75 * vy;
				vlooplist[tlist[triangle][j]][2] = 0.75 * vz;
			}
		}
	}
}

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


void myObjType::loopSubdivide(){
	std::map<std::pair<int, int>, int> odds =  std::map<std::pair<int, int>, int>();
	getNeighbours();
	int odd[3];
	int triangles = tcount;
	for (int i = 1; i <= triangles; ++i) {
		for (int j = 0; j < 3; ++j) {
			odd[j]=computeOdd(i,j,odds);
		}
		computeEven(i);
		mergeTriangles(i,odd);
	}
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
	orientTriangles();
	computeInfo();
}

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
	tlooplist[tcount][0] = even2;
	tlooplist[tcount][1] = odd2;
	tlooplist[tcount][2] = odd1;
	normal = computeTriangleNormal(vlooplist[even2], vlooplist[odd2], vlooplist[odd1]);
	nx = normal.at(0);
	ny = normal.at(1);
	nz = normal.at(2);
	if (nx * nlist[i][0] + ny * nlist[i][1] + nz * nlist[i][2] < 0) {
		tlooplist[tcount][0] = odd2;
		tlooplist[tcount][1] = even2;
		tlooplist[tcount][2] = odd1;
	}
	++tcount;
	tlooplist[tcount][0] = even1;
	tlooplist[tcount][1] = odd1;
	tlooplist[tcount][2] = odd3;
	normal = computeTriangleNormal(vlooplist[even1], vlooplist[odd1], vlooplist[odd3]);
	nx = normal.at(0);
	ny = normal.at(1);
	nz = normal.at(2);
	if (nx * nlist[i][0] + ny * nlist[i][1] + nz * nlist[i][2] < 0) {
		tlooplist[tcount][0] = odd1;
		tlooplist[tcount][1] = even1;
		tlooplist[tcount][2] = odd3;
	}
	++tcount;
	tlooplist[tcount][0] = even3;
	tlooplist[tcount][1] = odd3;
	tlooplist[tcount][2] = odd2;
	normal = computeTriangleNormal(vlooplist[even3], vlooplist[odd3], vlooplist[odd2]);
	nx = normal.at(0);
	ny = normal.at(1);
	nz = normal.at(2);
	if (nx * nlist[i][0] + ny * nlist[i][1] + nz * nlist[i][2] < 0) {
		tlooplist[tcount][0] = odd3;
		tlooplist[tcount][1] = even3;
		tlooplist[tcount][2] = odd2;
	}
}

void myObjType::load3DS(char* p_filename) {

}